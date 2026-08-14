"""M2 normalization pipeline consuming validated M1 artifacts."""

import json
import shutil
import subprocess
from collections import Counter
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from typing import Any
from uuid import uuid4

import polars as pl
from pydantic import BaseModel, ConfigDict

from genome_evidence import __version__
from genome_evidence.domain.variants import Variant
from genome_evidence.ingest.twenty_three_and_me import _validate_private_output
from genome_evidence.normalization.models import (
    CanonicalGenotype,
    LiftStatus,
    MappingCandidate,
    MappingOutcome,
    ObservationMapping,
    ReferenceValidation,
    StrandTransform,
)
from genome_evidence.normalization.resources import (
    FastaReferenceProvider,
    JsonLiftoverProvider,
    JsonMarkerProvider,
    canonical_assembly,
    canonical_chromosome,
)
from genome_evidence.qc.models import CallState, RawGenotypeObservation

ALGORITHM_VERSION = "m2-snv-1"


class NormalizationConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    marker_definitions: Path
    target_reference: Path
    target_build: str = "GRCh38"
    marker_version: str = "1"
    reference_version: str = "1"
    liftover: Path | None = None
    liftover_version: str = "1"
    source_build_override: str | None = None


class NormalizationResult(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    run_id: str
    output_directory: Path
    mappings: tuple[ObservationMapping, ...]
    variants: tuple[Variant, ...]
    genotypes: tuple[CanonicalGenotype, ...]
    candidates: tuple[MappingCandidate, ...]
    manifest: dict[str, Any]


def _hash(p: Path) -> str:
    return sha256(p.read_bytes()).hexdigest()


def _dump(v: Any) -> bytes:
    return (json.dumps(v, indent=2, sort_keys=True, default=str) + "\n").encode()


def _stable(prefix: str, v: Any) -> str:
    return (
        prefix
        + "-"
        + sha256(
            json.dumps(v, sort_keys=True, default=str, separators=(",", ":")).encode()
        ).hexdigest()
    )


def _transform(v: str, s: StrandTransform) -> str:
    c = v.translate(str.maketrans("ACGT", "TGCA"))
    return (
        c[::-1]
        if s == StrandTransform.REVERSE_COMPLEMENT
        else c
        if s == StrandTransform.COMPLEMENT
        else v
    )


def normalize_m1_run(m1: Path, output: Path, config: NormalizationConfig) -> NormalizationResult:
    required = [m1 / "manifest.json", m1 / "source_metadata.json", m1 / "observations.parquet"]
    if not all(p.is_file() for p in required):
        raise ValueError("incomplete M1 run directory")
    manifest = json.loads(required[0].read_text())
    meta = json.loads(required[1].read_text())
    for name, expected in manifest.get("artifacts", {}).items():
        if not (m1 / name).is_file() or _hash(m1 / name) != expected:
            raise ValueError(f"M1 artifact checksum mismatch: {name}")
    if manifest.get("run_id") != meta.get("run_id"):
        raise ValueError("incompatible M1 run identity")
    try:
        observations = tuple(
            RawGenotypeObservation.model_validate(x)
            for x in pl.read_parquet(required[2]).to_dicts()
        )
    except Exception as e:
        raise ValueError("incompatible M1 observations artifact") from e
    if any(x.ingestion_run_id != manifest["run_id"] for x in observations):
        raise ValueError("incompatible M1 observation identity")
    _validate_private_output(output)
    if output.exists():
        raise FileExistsError("normalization output already exists")
    marker = JsonMarkerProvider(config.marker_definitions, version=config.marker_version)
    target = canonical_assembly(config.target_build)
    if target is None:
        raise ValueError("unsupported target assembly")
    reference = FastaReferenceProvider(config.target_reference, target, config.reference_version)
    source_token = config.source_build_override or meta["resolved_build"]
    source = canonical_assembly(source_token)
    lift = (
        JsonLiftoverProvider(config.liftover, source, target, config.liftover_version)
        if source and source != target and config.liftover
        else None
    )
    confighash = sha256(_dump(config.model_dump(mode="json"))).hexdigest()
    run = str(uuid4())
    started = datetime.now(UTC)
    mappings = []
    variants = {}
    genotypes = []
    candidates = []
    for o in observations:
        obsref = _stable(
            "obs",
            [
                manifest["run_id"],
                manifest.get("input_sha256"),
                o.source_line_number,
                o.sample_id,
                o.source_marker_id,
                o.source_chromosome,
                o.source_position,
                o.raw_genotype,
            ],
        )
        common = dict(
            mapping_id=_stable("map", [obsref, confighash, ALGORITHM_VERSION]),
            observation_reference=obsref,
            normalization_run_id=run,
            source_assembly_token=source_token,
            resolved_source_assembly=source,
            source_chromosome=o.source_chromosome,
            source_position=o.source_position,
            target_assembly=target,
        )
        all_defs = marker.definitions(o.source_marker_id)
        defs = tuple(d for d in all_defs if canonical_assembly(d.assembly) == source)
        outcome = MappingOutcome.UNMAPPED
        reason: str | None = "MARKER_DEFINITION_ABSENT"
        strand = StrandTransform.UNKNOWN
        ls = LiftStatus.NOT_REQUIRED
        rv = ReferenceValidation.NOT_ATTEMPTED
        tc = None
        tp = None
        vid = None
        cids: tuple[str, ...] = ()
        correspondence: tuple[str, ...] = ()
        if source is None:
            reason = "MISSING_OR_UNKNOWN_BUILD"
        elif all_defs and not defs:
            outcome = MappingOutcome.UNSUPPORTED
            reason = "MARKER_DEFINITION_ASSEMBLY_MISMATCH"
        elif len(defs) != 1:
            if defs:
                outcome = MappingOutcome.AMBIGUOUS
                reason = "MARKER_DEFINITION_CONFLICT"
                made = [
                    MappingCandidate(
                        candidate_id=_stable("cand", [obsref, d.model_dump()]),
                        observation_reference=obsref,
                        normalization_run_id=run,
                        chromosome=canonical_chromosome(d.chromosome),
                        position=d.position,
                    )
                    for d in defs
                ]
                candidates += made
                cids = tuple(x.candidate_id for x in made)
        else:
            d = defs[0]
            strand = d.orientation
            ref = d.reference.upper()
            alt = d.alternate.upper()
            if any(not a or set(a) - set("ACGT") for a in (ref, alt)) or len(ref) != len(alt):
                outcome = MappingOutcome.UNSUPPORTED
                reason = "UNSUPPORTED_ALLELE_REPRESENTATION"
            elif not d.orientation_authoritative or strand in (
                StrandTransform.UNKNOWN,
                StrandTransform.AMBIGUOUS,
            ):
                outcome = MappingOutcome.AMBIGUOUS
                reason = (
                    "PALINDROMIC_AMBIGUITY"
                    if {ref, alt} in ({"A", "T"}, {"C", "G"})
                    else "STRAND_UNRESOLVED"
                )
                c = MappingCandidate(
                    candidate_id=_stable("cand", [obsref, d.model_dump()]),
                    observation_reference=obsref,
                    normalization_run_id=run,
                    chromosome=canonical_chromosome(d.chromosome),
                    position=d.position,
                )
                candidates.append(c)
                cids = (c.candidate_id,)
            else:
                tc = canonical_chromosome(d.chromosome)
                tp = d.position
                if source != target:
                    if lift is None:
                        reason = "NO_LIFTOVER_RESOURCE"
                        ls = LiftStatus.UNSUPPORTED
                    else:
                        lifted = lift.lift(tc, tp - 1)
                        made = [
                            MappingCandidate(
                                candidate_id=_stable("cand", [obsref, c, p]),
                                observation_reference=obsref,
                                normalization_run_id=run,
                                chromosome=c,
                                position=p + 1,
                            )
                            for c, p in lifted
                        ]
                        candidates += made
                        cids = tuple(x.candidate_id for x in made)
                        if not lifted:
                            reason = "NO_LIFTOVER_CANDIDATE"
                            ls = LiftStatus.FAILED
                        elif len(lifted) > 1:
                            outcome = MappingOutcome.AMBIGUOUS
                            reason = "MULTIPLE_LIFTOVER_CANDIDATES"
                            ls = LiftStatus.AMBIGUOUS
                        else:
                            tc, tp = lifted[0][0], lifted[0][1] + 1
                            ls = LiftStatus.SUCCESS
                if source == target or ls == LiftStatus.SUCCESS:
                    rv = (
                        ReferenceValidation.MATCH
                        if reference.sequence(tc, tp, len(ref)) == ref
                        else ReferenceValidation.MISMATCH
                    )
                    if rv == ReferenceValidation.MISMATCH:
                        outcome = MappingOutcome.FAILED
                        reason = "TARGET_REF_MISMATCH"
                    else:
                        v = Variant(
                            assembly=target,
                            chromosome=tc,
                            position=tp,
                            reference=ref,
                            alternate=alt,
                            rsid=d.rsid,
                        )
                        vid = _stable("var", [target, tc, tp, ref, alt])
                        variants[vid] = v
                        outcome = MappingOutcome.MAPPED
                        reason = None
                        if o.call_status == CallState.CALLED and len(o.raw_genotype) in (
                            len(ref),
                            2 * len(ref),
                        ):
                            alleles = tuple(
                                o.raw_genotype.upper()[i : i + len(ref)]
                                for i in range(0, len(o.raw_genotype), len(ref))
                            )
                            transformed = tuple(_transform(x, strand) for x in alleles)
                            correspondence = transformed
                            if all(x in (ref, alt) for x in transformed):
                                genotypes.append(
                                    CanonicalGenotype(
                                        genotype_id=_stable("gt", [obsref, vid, transformed]),
                                        observation_reference=obsref,
                                        normalization_run_id=run,
                                        variant_id=vid,
                                        alleles=transformed,
                                        ploidy=len(transformed),
                                        call_status=CallState.CALLED,
                                    )
                                )
        mappings.append(
            ObservationMapping(
                **common,
                outcome=outcome,
                reason=reason,
                target_chromosome=tc,
                target_position=tp,
                variant_id=vid,
                strand_transform=strand,
                liftover_required=source is not None and source != target,
                liftover_status=ls,
                reference_validation=rv,
                candidate_ids=cids,
                source_to_canonical_alleles=correspondence,
            )
        )
    temp = output.with_name(f".{output.name}.{run}.tmp")
    temp.mkdir(parents=True)
    try:
        tables = {
            "variants.parquet": [
                dict(variant_id=k, **v.model_dump()) for k, v in sorted(variants.items())
            ],
            "observation_mappings.parquet": [x.model_dump(mode="json") for x in mappings],
            "canonical_genotypes.parquet": [x.model_dump(mode="json") for x in genotypes],
            "mapping_candidates.parquet": [x.model_dump(mode="json") for x in candidates],
        }
        hashes = {}
        for name, rows in tables.items():
            pl.DataFrame(rows).write_parquet(temp / name) if rows else pl.DataFrame().write_parquet(
                temp / name
            )
            hashes[name] = _hash(temp / name)
        qc = {
            "outcomes": dict(Counter(x.outcome.value for x in mappings)),
            "reasons": dict(Counter(x.reason for x in mappings if x.reason)),
            "strand_transformations": dict(Counter(x.strand_transform.value for x in mappings)),
            "liftover": dict(Counter(x.liftover_status.value for x in mappings)),
            "reference_validation": dict(Counter(x.reference_validation.value for x in mappings)),
            "source_chromosomes": dict(Counter(x.source_chromosome for x in mappings)),
            "call_state": dict(Counter(x.call_status.value for x in observations)),
        }
        for name, data in [
            ("normalization_qc.json", _dump(qc)),
            (
                "normalization_report.md",
                (
                    "# Canonical normalization QC\n\n"
                    "> Aggregate results only; no genotypes or clinical interpretation.\n"
                    + "".join(f"- {k}: {v}\n" for k, v in sorted(qc["outcomes"].items()))
                ).encode(),
            ),
        ]:
            (temp / name).write_bytes(data)
            hashes[name] = _hash(temp / name)
        resources = [marker.identity.model_dump(), reference.identity.model_dump()] + (
            [lift.identity.model_dump()] if lift else []
        )
        md = {
            "run_id": run,
            "m1_run_id": manifest["run_id"],
            "m1_manifest_sha256": _hash(required[0]),
            "m1_observations_sha256": _hash(required[2]),
            "configuration": config.model_dump(mode="json"),
            "configuration_hash": confighash,
            "algorithm": ALGORITHM_VERSION,
            "package_version": __version__,
            "git_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], capture_output=True, text=True
            ).stdout.strip(),
            "source_assembly_token": source_token,
            "resolved_source_assembly": source,
            "target_assembly": target,
            "resources": resources,
            "observation_references": [x.observation_reference for x in mappings],
            "started_at": started.isoformat(),
            "completed_at": datetime.now(UTC).isoformat(),
        }
        (temp / "normalization_metadata.json").write_bytes(_dump(md))
        hashes["normalization_metadata.json"] = _hash(temp / "normalization_metadata.json")
        outmanifest = {
            "schema_version": 1,
            "run_id": run,
            "m1_run_id": manifest["run_id"],
            "artifacts": hashes,
        }
        (temp / "manifest.json").write_bytes(_dump(outmanifest))
        temp.rename(output)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    return NormalizationResult(
        run_id=run,
        output_directory=output,
        mappings=tuple(mappings),
        variants=tuple(variants.values()),
        genotypes=tuple(genotypes),
        candidates=tuple(candidates),
        manifest=outmanifest,
    )
