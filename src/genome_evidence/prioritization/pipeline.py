"""Deterministic, offline M4 manual-review routing."""

import json
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from typing import Any, cast
from uuid import uuid4

import polars as pl

from genome_evidence import __version__
from genome_evidence.evidence.models import (
    AssertionConditionRelationship,
    ExternalAssertion,
    ExternalCondition,
    ExternalSourceSnapshot,
    ExternalVariantRepresentation,
    LinkOutcome,
    VariantEvidenceLink,
)
from genome_evidence.ingest.twenty_three_and_me import _validate_private_output
from genome_evidence.normalization.models import CanonicalGenotype

from .models import (
    CandidateAssertionLink,
    CandidateEligibility,
    ClinicalPrioritizationConfig,
    ClinicalReviewCandidate,
    GenotypeEvidenceState,
    GenotypeRowEvidence,
    PolicyIdentity,
    PrioritizationExclusion,
    PrioritizationPolicy,
    PrioritizationResult,
    PriorityRationale,
    ReviewPriorityBand,
    SourceReviewLevel,
    SourceTermRoute,
    VariantEvidenceProfile,
)

ALGORITHM_VERSION = "m4-clinvar-germline-review-routing-1"
UNRESOLVED = (
    "phenotype_concordance_not_assessed",
    "family_history_not_assessed",
    "inheritance_fit_not_assessed",
    "penetrance_not_assessed",
    "actionability_not_assessed",
    "clinical_confirmation_not_assessed",
    "family_segregation_not_assessed",
)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, default=str) + "\n").encode()


def _stable(prefix: str, value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return f"{prefix}-{sha256(raw).hexdigest()}"


def _manifest(directory: Path, label: str) -> dict[str, Any]:
    path = directory / "manifest.json"
    if not path.is_file():
        raise ValueError(f"missing {label} manifest")
    try:
        result = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as error:
        raise ValueError(f"invalid {label} manifest") from error
    if result.get("schema_version") != 1 or not result.get("run_id"):
        raise ValueError(f"unsupported {label} manifest schema")
    if not isinstance(result.get("artifacts"), dict):
        raise ValueError(f"invalid {label} artifact registry")
    for name, expected in result["artifacts"].items():
        artifact = directory / name
        if not artifact.is_file() or _hash(artifact) != expected:
            raise ValueError(f"{label} artifact checksum mismatch: {name}")
    return cast(dict[str, Any], result)


def _rows(path: Path, model: type[Any], label: str) -> list[Any]:
    try:
        return [model.model_validate(row) for row in pl.read_parquet(path).to_dicts()]
    except Exception as error:
        raise ValueError(f"incompatible {label} Parquet schema") from error


def _unique(rows: list[Any], attr: str, label: str) -> None:
    values = [getattr(row, attr) for row in rows]
    if len(values) != len(set(values)):
        raise ValueError(f"duplicate {label} identifiers")


def _policy(config: ClinicalPrioritizationConfig) -> tuple[PrioritizationPolicy, PolicyIdentity]:
    raw = config.policy_path.read_bytes()
    try:
        parsed = json.loads(raw)
        policy = PrioritizationPolicy.model_validate(parsed)
    except Exception as error:
        raise ValueError("invalid prioritization policy") from error
    if policy.supported_analysis_context != config.analysis_context:
        raise ValueError("policy and analysis context mismatch")
    canonical = policy.model_dump(mode="json")
    return policy, PolicyIdentity(
        policy_id=policy.policy_id,
        version=policy.version,
        file_sha256=sha256(raw).hexdigest(),
        file_size_bytes=len(raw),
        configuration_hash=sha256(_json(canonical)).hexdigest(),
        parsed_configuration=canonical,
    )


def _normalize(value: str) -> str:
    return " ".join(value.casefold().split())


def _genotype_state(
    genotypes: list[CanonicalGenotype], reference: str, alternate: str
) -> tuple[GenotypeEvidenceState, tuple[GenotypeRowEvidence, ...]]:
    rows = []
    multisets = set()
    for genotype in sorted(genotypes, key=lambda x: x.genotype_id):
        if any(allele not in {reference, alternate} for allele in genotype.alleles):
            raise ValueError("canonical genotype contains an allele outside its REF/ALT set")
        multisets.add(tuple(sorted(genotype.alleles)))
        rows.append(
            GenotypeRowEvidence(
                genotype_id=genotype.genotype_id,
                observation_reference=genotype.observation_reference,
                ploidy=genotype.ploidy,
                alt_allele_count=sum(x == alternate for x in genotype.alleles),
            )
        )
    if not rows:
        state = GenotypeEvidenceState.NO_CANONICAL_CALLED_GENOTYPE
    elif len(multisets) > 1:
        state = (
            GenotypeEvidenceState.DISCORDANT_CALLED_ROWS_WITH_ALT
            if any(row.alt_allele_count for row in rows)
            else GenotypeEvidenceState.OBSERVED_REFERENCE_ONLY
        )
    elif not any(row.alt_allele_count for row in rows):
        state = GenotypeEvidenceState.OBSERVED_REFERENCE_ONLY
    elif len(rows) == 1:
        state = GenotypeEvidenceState.SINGLE_OBSERVED_ALT_CALL
    else:
        state = GenotypeEvidenceState.CONCORDANT_OBSERVED_ALT_CALLS
    return state, tuple(rows)


def _write(path: Path, records: list[Any], schema: dict[str, Any]) -> None:
    rows = [record.model_dump(mode="json") for record in records]
    (pl.DataFrame(rows) if rows else pl.DataFrame(schema=schema)).write_parquet(path)


def prioritize_clinical_variants(
    normalization_directory: Path,
    evidence_directory: Path,
    annotation_directory: Path,
    output_directory: Path,
    config: ClinicalPrioritizationConfig,
) -> PrioritizationResult:
    """Build an auditable review queue; never classify a variant or subject."""
    directories = [normalization_directory, evidence_directory, annotation_directory]
    normalization_directory, evidence_directory, annotation_directory = [
        path.resolve() for path in directories
    ]
    output_directory = output_directory.resolve()
    m2 = _manifest(normalization_directory, "M2")
    evidence = _manifest(evidence_directory, "evidence")
    annotation = _manifest(annotation_directory, "annotation")
    if (
        annotation.get("m2_run_id") != m2["run_id"]
        or annotation.get("evidence_run_id") != evidence["run_id"]
    ):
        raise ValueError("cross-run M2/evidence/annotation identity mismatch")
    policy, identity = _policy(config)
    metadata = json.loads((evidence_directory / "external_source_metadata.json").read_text())
    snapshot = ExternalSourceSnapshot.model_validate(metadata["snapshot"])
    if snapshot.source_namespace != policy.source_namespace or snapshot.dataset != policy.dataset:
        raise ValueError("policy does not support supplied source snapshot")

    variant_rows = pl.read_parquet(normalization_directory / "variants.parquet").to_dicts()
    if not {"variant_id", "assembly", "chromosome", "position", "reference", "alternate"} <= set(
        pl.read_parquet(normalization_directory / "variants.parquet").columns
    ):
        raise ValueError("incompatible M2 variants Parquet schema")
    variant_ids = [row["variant_id"] for row in variant_rows]
    if len(variant_ids) != len(set(variant_ids)):
        raise ValueError("duplicate M2 variant identifiers")
    exact_keys = [
        (r["assembly"], r["chromosome"], r["position"], r["reference"], r["alternate"])
        for r in variant_rows
    ]
    if len(exact_keys) != len(set(exact_keys)):
        raise ValueError("duplicate exact M2 allele keys are ambiguous")
    genotypes = _rows(
        normalization_directory / "canonical_genotypes.parquet", CanonicalGenotype, "M2 genotypes"
    )
    _unique(genotypes, "genotype_id", "genotype")
    if any(g.normalization_run_id != m2["run_id"] for g in genotypes):
        raise ValueError("M2 genotype run identity mismatch")
    if any(g.variant_id not in set(variant_ids) for g in genotypes):
        raise ValueError("M2 genotype references an unknown variant")
    representations = _rows(
        evidence_directory / "external_variant_representations.parquet",
        ExternalVariantRepresentation,
        "representations",
    )
    assertions = _rows(
        evidence_directory / "external_assertions.parquet", ExternalAssertion, "assertions"
    )
    conditions = _rows(
        evidence_directory / "external_conditions.parquet", ExternalCondition, "conditions"
    )
    condition_links = _rows(
        evidence_directory / "assertion_conditions.parquet",
        AssertionConditionRelationship,
        "assertion conditions",
    )
    links = _rows(
        annotation_directory / "variant_evidence_links.parquet",
        VariantEvidenceLink,
        "evidence links",
    )
    for rows, attr, label in (
        (representations, "representation_id", "representation"),
        (assertions, "assertion_instance_id", "assertion"),
        (conditions, "condition_id", "condition"),
        (links, "link_id", "link"),
    ):
        _unique(rows, attr, label)
    rep_ids, assertion_ids, condition_ids = (
        {x.representation_id for x in representations},
        {x.assertion_instance_id for x in assertions},
        {x.condition_id for x in conditions},
    )
    if any(a.representation_id not in rep_ids for a in assertions) or any(
        link.representation_id not in rep_ids for link in links
    ):
        raise ValueError("evidence record references an unknown representation")
    if any(
        link.outcome == LinkOutcome.MATCHED and link.variant_id not in set(variant_ids)
        for link in links
    ):
        raise ValueError("matched link references an unknown M2 variant")
    if any(
        x.assertion_instance_id not in assertion_ids or x.condition_id not in condition_ids
        for x in condition_links
    ):
        raise ValueError("assertion-condition relationship has an unknown endpoint")

    route_index = {
        _normalize(term): route for route, terms in policy.term_routes.items() for term in terms
    }
    by_variant: dict[str, list[VariantEvidenceLink]] = defaultdict(list)
    for link in links:
        if link.outcome == LinkOutcome.MATCHED and link.variant_id:
            by_variant[link.variant_id].append(link)
    assertion_by_rep: dict[str, list[ExternalAssertion]] = defaultdict(list)
    for assertion in assertions:
        assertion_by_rep[assertion.representation_id].append(assertion)
    condition_by_assertion: dict[str, list[str]] = defaultdict(list)
    for relation in condition_links:
        condition_by_assertion[relation.assertion_instance_id].append(relation.condition_id)
    condition_map = {x.condition_id: x for x in conditions}
    genotype_by_variant: dict[str, list[CanonicalGenotype]] = defaultdict(list)
    for genotype in genotypes:
        genotype_by_variant[genotype.variant_id].append(genotype)

    profiles: list[VariantEvidenceProfile] = []
    candidates: list[ClinicalReviewCandidate] = []
    candidate_links: list[CandidateAssertionLink] = []
    rationales: list[PriorityRationale] = []
    exclusions: list[PrioritizationExclusion] = []
    for variant in sorted(variant_rows, key=lambda x: x["variant_id"]):
        matched = by_variant.get(variant["variant_id"], [])
        if not matched:
            continue
        reps = sorted({x.representation_id for x in matched})
        relevant = sorted(
            (a for rep in reps for a in assertion_by_rep[rep]),
            key=lambda x: x.assertion_instance_id,
        )
        state, genotype_rows = _genotype_state(
            genotype_by_variant[variant["variant_id"]], variant["reference"], variant["alternate"]
        )
        routed = [
            (a, term, route_index.get(_normalize(term), SourceTermRoute.UNMAPPED))
            for a in relevant
            for term in a.source_classification_terms
        ]
        active = [
            (a, t, r)
            for a, t, r in routed
            if policy.record_status_behavior.get(a.source_record_status.value, "inactive")
            == "active"
        ]
        active_germline = [
            (a, t, r) for a, t, r in active if a.classification_type.value == "germline"
        ]
        source_conflict = any(
            "conflict" in _normalize(term)
            for a, term, _ in routed
            if a.assertion_level.value == "source_computed_aggregate"
            and a.classification_type.value == "germline"
        )
        profile_id = _stable(
            "profile",
            [
                m2["run_id"],
                variant["variant_id"],
                reps,
                identity.configuration_hash,
                ALGORITHM_VERSION,
            ],
        )
        cids = sorted(
            {cid for a in relevant for cid in condition_by_assertion[a.assertion_instance_id]}
        )
        profile = VariantEvidenceProfile(
            profile_id=profile_id,
            m2_run_id=m2["run_id"],
            variant_id=variant["variant_id"],
            assembly=variant["assembly"],
            chromosome=variant["chromosome"],
            position=variant["position"],
            reference=variant["reference"],
            alternate=variant["alternate"],
            genotype_state=state,
            genotype_rows=genotype_rows,
            representation_ids=tuple(reps),
            assertion_ids=tuple(a.assertion_instance_id for a in relevant),
            scv_assertion_ids=tuple(a.assertion_instance_id for a in relevant if a.scv_accession),
            vcv_assertion_ids=tuple(
                a.assertion_instance_id for a in relevant if not a.scv_accession
            ),
            assertion_levels=tuple(a.assertion_level.value for a in relevant),
            classification_types=tuple(a.classification_type.value for a in relevant),
            source_terms=tuple(term for _, term, _ in routed),
            source_term_routes=tuple(route for _, _, route in routed),
            source_review_statuses=tuple(a.source_review_status or "unknown" for a in relevant),
            source_review_levels=tuple(
                policy.source_review_status_mappings.get(
                    _normalize(a.source_review_status or ""), SourceReviewLevel.UNKNOWN
                )
                for a in relevant
            ),
            source_record_statuses=tuple(a.source_record_status.value for a in relevant),
            submitters=tuple(a.submitter_name for a in relevant if a.submitter_name),
            condition_ids=tuple(cids),
            condition_names=tuple(condition_map[c].source_name for c in cids),
            date_last_evaluated=tuple(
                str(a.date_last_evaluated) if a.date_last_evaluated else "unknown" for a in relevant
            ),
            source_snapshot_id=snapshot.snapshot_id,
            source_release_date=str(snapshot.release_date),
            source_reported_conflict=source_conflict,
            submission_term_diversity=tuple(
                sorted(
                    {
                        t
                        for a, t, _ in routed
                        if a.scv_accession and a.classification_type.value == "germline"
                    }
                )
            ),
            missing_indicators=tuple(
                sorted({"date_last_evaluated" for a in relevant if a.date_last_evaluated is None})
            ),
            unresolved_assessments=UNRESOLVED,
        )
        profiles.append(profile)
        alt_eligible = state in {
            GenotypeEvidenceState.SINGLE_OBSERVED_ALT_CALL,
            GenotypeEvidenceState.CONCORDANT_OBSERVED_ALT_CALLS,
            GenotypeEvidenceState.DISCORDANT_CALLED_ROWS_WITH_ALT,
        }
        if not alt_eligible:
            band, eligibility, reason = (
                ReviewPriorityBand.NOT_ELIGIBLE,
                CandidateEligibility.NOT_ELIGIBLE,
                state.value,
            )
        elif any(r == SourceTermRoute.HIGH_ATTENTION for _, _, r in active_germline):
            band, eligibility, reason = (
                ReviewPriorityBand.REVIEW_FIRST,
                CandidateEligibility.DATA_CONFLICT
                if state == GenotypeEvidenceState.DISCORDANT_CALLED_ROWS_WITH_ALT
                else CandidateEligibility.ELIGIBLE,
                "active_high_attention_source_term",
            )
        elif source_conflict or any(
            r
            in {
                SourceTermRoute.RISK_CONTEXT,
                SourceTermRoute.UNCERTAIN,
                SourceTermRoute.OTHER_CONTEXT,
                SourceTermRoute.UNMAPPED,
            }
            for _, _, r in active_germline
        ):
            band, eligibility, reason = (
                ReviewPriorityBand.REVIEW_NEXT,
                CandidateEligibility.DATA_CONFLICT
                if state == GenotypeEvidenceState.DISCORDANT_CALLED_ROWS_WITH_ALT
                else CandidateEligibility.ELIGIBLE,
                "active_germline_review_context",
            )
        elif active_germline and all(
            r == SourceTermRoute.BENIGN_LIKE for _, _, r in active_germline
        ):
            band, eligibility, reason = (
                ReviewPriorityBand.NOT_PRIORITIZED,
                CandidateEligibility.ELIGIBLE,
                "benign_like_only",
            )
        else:
            band, eligibility, reason = (
                ReviewPriorityBand.CONTEXT_ONLY,
                CandidateEligibility.ELIGIBLE,
                "no_active_germline_assertion",
            )
        candidate_id = _stable("candidate", [profile_id, band.value])
        candidate = ClinicalReviewCandidate(
            candidate_id=candidate_id,
            profile_id=profile_id,
            variant_id=variant["variant_id"],
            priority_band=band,
            eligibility=eligibility,
            ordering_components=(band.value, "source_review_status_contextual", profile_id),
        )
        candidates.append(candidate)
        for assertion, term, route in routed:
            candidate_links.append(
                CandidateAssertionLink(
                    link_id=_stable(
                        "candidate-assertion", [candidate_id, assertion.assertion_instance_id, term]
                    ),
                    candidate_id=candidate_id,
                    assertion_instance_id=assertion.assertion_instance_id,
                    assertion_level=assertion.assertion_level.value,
                    classification_type=assertion.classification_type.value,
                    source_term=term,
                    source_term_route=route,
                    active_for_routing=(assertion, term, route) in active_germline,
                )
            )
        rationale_type = (
            "excludes"
            if band == ReviewPriorityBand.NOT_ELIGIBLE
            else "promotes"
            if band in {ReviewPriorityBand.REVIEW_FIRST, ReviewPriorityBand.REVIEW_NEXT}
            else "contextualizes"
        )
        rationales.append(
            PriorityRationale(
                rationale_id=_stable("rationale", [candidate_id, reason]),
                candidate_id=candidate_id,
                profile_id=profile_id,
                policy_rule_id=reason,
                policy_rule_version=policy.version,
                rationale_type=rationale_type,
                reason_code=reason,
                assertion_ids=tuple(a.assertion_instance_id for a, _, _ in active_germline),
                genotype_ids=tuple(x.genotype_id for x in genotype_rows),
                observation_references=tuple(x.observation_reference for x in genotype_rows),
                explanation=(
                    f"Policy {policy.policy_id} routed this source-linked record to "
                    f"{band.value}; the band is only review ordering."
                ),
            )
        )
        if eligibility == CandidateEligibility.DATA_CONFLICT:
            rationales.append(
                PriorityRationale(
                    rationale_id=_stable("rationale", [candidate_id, "data_conflict"]),
                    candidate_id=candidate_id,
                    profile_id=profile_id,
                    policy_rule_id="discordant-called-rows",
                    policy_rule_version=policy.version,
                    rationale_type="warns",
                    reason_code="data_conflict_alt_carriage_unresolved",
                    genotype_ids=tuple(x.genotype_id for x in genotype_rows),
                    observation_references=tuple(x.observation_reference for x in genotype_rows),
                    explanation=(
                        "Called rows disagree; ALT carriage remains unresolved and no row was "
                        "selected by vote."
                    ),
                )
            )
        if band == ReviewPriorityBand.NOT_ELIGIBLE:
            exclusions.append(
                PrioritizationExclusion(
                    exclusion_id=_stable("exclusion", [profile_id, reason]),
                    profile_id=profile_id,
                    candidate_id=candidate_id,
                    reason_code=reason,
                )
            )

    band_order = {
        b: i
        for i, b in enumerate(
            (
                ReviewPriorityBand.REVIEW_FIRST,
                ReviewPriorityBand.REVIEW_NEXT,
                ReviewPriorityBand.CONTEXT_ONLY,
                ReviewPriorityBand.NOT_PRIORITIZED,
                ReviewPriorityBand.NOT_ELIGIBLE,
            )
        )
    }
    candidates.sort(key=lambda x: (band_order[x.priority_band], x.profile_id))
    run_id, started = str(uuid4()), datetime.now(UTC)
    _validate_private_output(output_directory)
    if output_directory.exists():
        raise FileExistsError("prioritization output already exists")
    temporary = output_directory.with_name(f".{output_directory.name}.{run_id}.tmp")
    temporary.mkdir(parents=True)
    try:
        specs: dict[str, tuple[list[Any], dict[str, Any]]] = {
            "variant_evidence_profiles.parquet": (profiles, {"profile_id": pl.String}),
            "prioritization_candidates.parquet": (candidates, {"candidate_id": pl.String}),
            "candidate_assertion_links.parquet": (candidate_links, {"link_id": pl.String}),
            "priority_rationales.parquet": (rationales, {"rationale_id": pl.String}),
            "prioritization_exclusions.parquet": (exclusions, {"exclusion_id": pl.String}),
        }
        artifacts = {}
        for name, (records, schema) in specs.items():
            _write(temporary / name, records, schema)
            artifacts[name] = _hash(temporary / name)
        configuration = config.model_dump(mode="json")
        meta = {
            "run_id": run_id,
            "m2_run_id": m2["run_id"],
            "evidence_run_id": evidence["run_id"],
            "annotation_run_id": annotation["run_id"],
            "upstream_manifest_hashes": {
                "m2": _hash(normalization_directory / "manifest.json"),
                "evidence": _hash(evidence_directory / "manifest.json"),
                "annotation": _hash(annotation_directory / "manifest.json"),
            },
            "policy_identity": identity.model_dump(mode="json"),
            "configuration": configuration,
            "configuration_hash": sha256(_json(configuration)).hexdigest(),
            "algorithm_version": ALGORITHM_VERSION,
            "package_version": __version__,
            "git_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False
            ).stdout.strip(),
            "started_at": started.isoformat(),
            "completed_at": datetime.now(UTC).isoformat(),
            "analysis_context": config.analysis_context.value,
            "input_counts": {
                "variants": len(variant_rows),
                "exact_linked_profiles": len(profiles),
                "assertions": len(assertions),
            },
        }
        qc = {
            "priority_bands": dict(Counter(x.priority_band.value for x in candidates)),
            "genotype_states": dict(Counter(x.genotype_state.value for x in profiles)),
            "data_conflicts": sum(
                x.eligibility == CandidateEligibility.DATA_CONFLICT for x in candidates
            ),
        }
        for name, value in (("prioritization_metadata.json", meta), ("prioritization_qc.json", qc)):
            (temporary / name).write_bytes(_json(value))
            artifacts[name] = _hash(temporary / name)
        report = _report(candidates, profiles, candidate_links, rationales)
        (temporary / "prioritization_report.md").write_text(report)
        artifacts["prioritization_report.md"] = _hash(temporary / "prioritization_report.md")
        manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "m2_run_id": m2["run_id"],
            "evidence_run_id": evidence["run_id"],
            "annotation_run_id": annotation["run_id"],
            "policy_identity": identity.model_dump(mode="json"),
            "artifacts": artifacts,
        }
        (temporary / "manifest.json").write_bytes(_json(manifest))
        temporary.rename(output_directory)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return PrioritizationResult(
        run_id=run_id,
        output_directory=output_directory,
        policy_identity=identity,
        profiles=tuple(profiles),
        candidates=tuple(candidates),
        rationales=tuple(rationales),
        manifest=manifest,
    )


def _report(
    candidates: list[ClinicalReviewCandidate],
    profiles: list[VariantEvidenceProfile],
    links: list[CandidateAssertionLink],
    rationales: list[PriorityRationale],
) -> str:
    lines = [
        "# Research clinical-variant review queue",
        "",
        "> **Research review queue only. Source classifications are not project conclusions.**",
        "> Observed DTC-array calls are not clinically confirmed. A review band is not a "
        "disease-risk, pathogenicity, urgency, confidence, or actionability score.",
        "> Phenotype, family history, inheritance, penetrance, intervention evidence, "
        "clinical confirmation, and family segregation were not assessed.",
        "> Absence from this queue is not a negative genetic test; array content is not "
        "comprehensive genome coverage.",
        "> No medical decisions should be based on this report alone.",
        "",
    ]
    by_profile = {p.profile_id: p for p in profiles}
    for candidate in candidates:
        p = by_profile[candidate.profile_id]
        cl = [x for x in links if x.candidate_id == candidate.candidate_id]
        rr = [x.reason_code for x in rationales if x.candidate_id == candidate.candidate_id]
        lines += [
            f"## {p.assembly}:{p.chromosome}:{p.position}:{p.reference}>{p.alternate}",
            f"- Review band: `{candidate.priority_band.value}` (manual routing only)",
            f"- Genotype evidence state: `{p.genotype_state.value}`",
            "- ClinVar records: "
            + (", ".join(sorted(set(p.scv_assertion_ids + p.vcv_assertion_ids))) or "none"),
            f"- Source terms: {', '.join(x.source_term for x in cl) or 'none'}",
            f"- Review statuses: {', '.join(p.source_review_statuses)}",
            f"- Assertion levels: {', '.join(p.assertion_levels)}",
            f"- Conditions: {', '.join(p.condition_names) or 'not provided'}",
            f"- Source-reported conflict: {str(p.source_reported_conflict).lower()}",
            f"- Evaluation dates: {', '.join(p.date_last_evaluated)}",
            f"- Policy reasons: {', '.join(rr)}",
            f"- Unresolved: {', '.join(p.unresolved_assessments)}",
            "",
        ]
    return "\n".join(lines) + "\n"
