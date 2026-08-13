"""Offline ClinVar VCV ingestion and exact M2 allele linking."""

import gzip
import json
import shutil
import subprocess
from collections import Counter
from collections.abc import Iterator
from datetime import UTC, date, datetime
from hashlib import sha256
from pathlib import Path
from typing import IO, Any, cast
from uuid import uuid4
from xml.etree import ElementTree as ET

import polars as pl

from genome_evidence import __version__
from genome_evidence.ingest.twenty_three_and_me import _validate_private_output

from .models import (
    AnnotationConfig,
    AnnotationResult,
    AssertionConditionRelationship,
    AssertionLevel,
    AssertionRelationship,
    ClassificationType,
    ClinVarIngestionConfig,
    EvidenceIngestionResult,
    ExternalAssertion,
    ExternalCondition,
    ExternalSourceSnapshot,
    ExternalVariantRepresentation,
    LinkOutcome,
    LinkReason,
    SourceRecordStatus,
    VariantEvidenceLink,
)

PARSER_VERSION = "clinvar-vcv-xml-1"
LINKER_VERSION = "m3-exact-allele-1"


def _json(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, default=str) + "\n").encode()


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _stable(prefix: str, value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return f"{prefix}-{sha256(encoded).hexdigest()}"


def _local(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def _children(element: ET.Element, *names: str) -> Iterator[ET.Element]:
    wanted = set(names)
    return (node for node in element.iter() if _local(node) in wanted)


def _first_text(element: ET.Element, *names: str) -> str | None:
    for node in _children(element, *names):
        if node.text and node.text.strip():
            return node.text.strip()
    return None


def _parse_date(value: str | None) -> date | None:
    if not value:
        return None
    try:
        return date.fromisoformat(value[:10])
    except ValueError:
        return None


def _accession(element: ET.Element, fallback: str) -> tuple[str, int]:
    node = next(_children(element, "ClinVarAccession"), None)
    accession = (
        element.attrib.get("Accession")
        or (node.attrib.get("Acc") if node is not None else None)
        or fallback
    )
    raw_version = (
        element.attrib.get("Version")
        or (node.attrib.get("Version") if node is not None else None)
        or "1"
    )
    try:
        return accession, int(raw_version)
    except ValueError as error:
        raise ValueError(f"invalid accession version for {accession}") from error


def _classification_type(tag: str) -> ClassificationType:
    return {
        "GermlineClassification": ClassificationType.GERMLINE,
        "ClinicalSignificance": ClassificationType.GERMLINE,
        "SomaticClinicalImpact": ClassificationType.SOMATIC_CLINICAL_IMPACT,
        "OncogenicityClassification": ClassificationType.ONCOGENICITY,
    }.get(tag, ClassificationType.UNKNOWN)


def _classifications(element: ET.Element) -> list[tuple[ClassificationType, ET.Element]]:
    names = {
        "GermlineClassification",
        "ClinicalSignificance",
        "SomaticClinicalImpact",
        "OncogenicityClassification",
    }
    found = [(_classification_type(_local(n)), n) for n in element.iter() if _local(n) in names]
    # ClinicalAssertion wrappers can carry a direct Description in older VCV XML.
    return found or [(ClassificationType.UNKNOWN, element)]


def _terms(element: ET.Element) -> tuple[str, ...]:
    values = []
    for node in _children(element, "Description", "Classification", "ClinicalSignificance"):
        if node is element:
            continue
        if node.text and node.text.strip() and not list(node):
            values.append(node.text.strip())
    return tuple(dict.fromkeys(values))


def _record_status(term: str | None) -> SourceRecordStatus:
    normalized = (term or "").lower()
    if "current" in normalized:
        return SourceRecordStatus.CURRENT
    if "replac" in normalized:
        return SourceRecordStatus.REPLACED
    if "delet" in normalized:
        return SourceRecordStatus.DELETED
    return SourceRecordStatus.UNKNOWN


def _write_table(path: Path, rows: list[dict[str, Any]], schema: dict[str, Any]) -> None:
    frame = pl.DataFrame(rows) if rows else pl.DataFrame(schema=schema)
    frame.write_parquet(path)


def _prepare_output(output: Path, run_id: str) -> Path:
    _validate_private_output(output)
    if output.exists():
        raise FileExistsError("output directory already exists")
    temporary = output.with_name(f".{output.name}.{run_id}.tmp")
    temporary.mkdir(parents=True)
    return temporary


def _git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False
    ).stdout.strip()


def _open_xml(path: Path) -> IO[bytes]:
    # Ownership is transferred to the caller, which always enters this stream as a context manager.
    stream = gzip.open(path, "rb") if path.suffix.lower() == ".gz" else path.open("rb")  # noqa: SIM115
    return cast(IO[bytes], stream)


def ingest_clinvar_vcv(
    input_path: Path,
    output_directory: Path,
    config: ClinVarIngestionConfig | None = None,
) -> EvidenceIngestionResult:
    """Stream one explicit local ClinVar VCV XML snapshot into immutable artifacts."""
    config = config or ClinVarIngestionConfig()
    input_path, output_directory = input_path.resolve(), output_directory.resolve()
    if not input_path.is_file():
        raise FileNotFoundError("ClinVar source file does not exist")
    exact_hash, size = _hash(input_path), input_path.stat().st_size
    retrieved = config.retrieval_timestamp or datetime.now(UTC)
    run_id = str(uuid4())
    temporary = _prepare_output(output_directory, run_id)
    variants: list[ExternalVariantRepresentation] = []
    assertions: list[ExternalAssertion] = []
    relationships: list[AssertionRelationship] = []
    conditions: dict[str, ExternalCondition] = {}
    condition_links: list[AssertionConditionRelationship] = []
    root_attrs: dict[str, str] = {}
    release_date: date | None = None
    release_identity: str | None = None
    schema_identity: str | None = None
    try:
        try:
            with _open_xml(input_path) as stream:
                context = ET.iterparse(stream, events=("start", "end"))
                _, root = next(context)
                root_name = _local(root)
                if root_name not in {"ReleaseSet", "ClinVarVariationRelease"}:
                    raise ValueError(f"unsupported ClinVar VCV XML root: {root_name}")
                root_attrs = dict(root.attrib)
                raw_date = next(
                    (
                        root.attrib.get(k)
                        for k in ("Dated", "ReleaseDate", "releaseDate")
                        if root.attrib.get(k)
                    ),
                    None,
                )
                release_date = _parse_date(raw_date)
                release_identity = root.attrib.get("ReleaseID") or root.attrib.get("Release")
                schema_identity = next(
                    (v for k, v in root.attrib.items() if k.endswith("schemaLocation")), None
                ) or root.attrib.get("SchemaVersion")
                if (
                    release_identity
                    and config.release_identity_override
                    and release_identity != config.release_identity_override
                ):
                    raise ValueError("release identity override contradicts XML metadata")
                release_identity = release_identity or config.release_identity_override or raw_date
                if not release_identity or release_date is None:
                    raise ValueError(
                        "release identity/date absent; provide a non-contradictory release override"
                    )
                snapshot_id = _stable(
                    "snapshot", ["ClinVar", "VCV_XML", release_identity, exact_hash]
                )
                for event, archive in context:
                    if event != "end" or _local(archive) != "VariationArchive":
                        continue
                    vcv, vcv_version = _accession(archive, "")
                    if not vcv.startswith("VCV"):
                        raise ValueError("VariationArchive lacks a VCV accession")
                    classified = next(
                        _children(archive, "ClassifiedRecord", "IncludedRecord"), archive
                    )
                    allele = next(
                        _children(classified, "SimpleAllele", "Haplotype", "Genotype"), None
                    )
                    variation_type = _local(allele) if allele is not None else "unknown"
                    variation_id = archive.attrib.get("VariationID") or (
                        allele.attrib.get("VariationID") if allele is not None else None
                    )
                    allele_id = allele.attrib.get("AlleleID") if allele is not None else None
                    location_parent = allele if allele is not None else classified
                    location = next(_children(location_parent, "SequenceLocation"), None)
                    assembly = location.attrib.get("Assembly") if location is not None else None
                    chromosome = location.attrib.get("Chr") if location is not None else None
                    pos = None
                    if location is not None:
                        raw_pos = location.attrib.get("positionVCF") or location.attrib.get("start")
                        pos = int(raw_pos) if raw_pos and raw_pos.isdigit() else None
                    reference = (
                        location.attrib.get("referenceAlleleVCF") if location is not None else None
                    )
                    alternate = (
                        location.attrib.get("alternateAlleleVCF") if location is not None else None
                    )
                    representation_id = _stable(
                        "extvar",
                        [
                            snapshot_id,
                            vcv,
                            vcv_version,
                            assembly,
                            chromosome,
                            pos,
                            reference,
                            alternate,
                            variation_type,
                        ],
                    )
                    variants.append(
                        ExternalVariantRepresentation(
                            representation_id=representation_id,
                            snapshot_id=snapshot_id,
                            vcv_accession=vcv,
                            vcv_version=vcv_version,
                            variation_id=int(variation_id)
                            if variation_id and variation_id.isdigit()
                            else None,
                            allele_id=int(allele_id) if allele_id and allele_id.isdigit() else None,
                            variation_type=variation_type,
                            assembly=assembly,
                            chromosome=chromosome,
                            position=pos,
                            reference=reference,
                            alternate=alternate,
                            source_rsid=_first_text(classified, "XRef") if False else None,
                        )
                    )
                    made_aggregates: list[ExternalAssertion] = []
                    classification_container = next(
                        _children(classified, "Classifications"), classified
                    )
                    for ctype, node in _classifications(classification_container):
                        # Do not mistake classifications inside SCV submissions for aggregates.
                        if any(
                            node is nested
                            for ca in _children(classified, "ClinicalAssertion")
                            for nested in ca.iter()
                        ):
                            continue
                        made_aggregates.append(
                            _make_assertion(
                                node,
                                archive,
                                snapshot_id,
                                representation_id,
                                vcv,
                                vcv_version,
                                AssertionLevel.AGGREGATE,
                                ctype,
                                None,
                            )
                        )
                    if not made_aggregates:
                        made_aggregates.append(
                            _make_assertion(
                                classification_container,
                                archive,
                                snapshot_id,
                                representation_id,
                                vcv,
                                vcv_version,
                                AssertionLevel.AGGREGATE,
                                ClassificationType.UNKNOWN,
                                None,
                            )
                        )
                    assertions.extend(made_aggregates)
                    for clinical in _children(classified, "ClinicalAssertion"):
                        scv, scv_version = _accession(clinical, "")
                        if not scv.startswith("SCV"):
                            continue
                        submitter = next(_children(clinical, "Submitter"), None)
                        made: list[ExternalAssertion] = []
                        for ctype, node in _classifications(clinical):
                            item = _make_assertion(
                                node,
                                clinical,
                                snapshot_id,
                                representation_id,
                                scv,
                                scv_version,
                                AssertionLevel.SUBMITTED,
                                ctype,
                                submitter,
                                vcv_accession=vcv,
                            )
                            made.append(item)
                            for aggregate in made_aggregates:
                                relationships.append(
                                    AssertionRelationship(
                                        relationship_id=_stable(
                                            "arel",
                                            [
                                                item.assertion_instance_id,
                                                aggregate.assertion_instance_id,
                                            ],
                                        ),
                                        snapshot_id=snapshot_id,
                                        subject_assertion_id=item.assertion_instance_id,
                                        predicate="contributes_to_aggregate",
                                        object_assertion_id=aggregate.assertion_instance_id,
                                    )
                                )
                            for trait in _children(clinical, "Trait", "Condition"):
                                name = _first_text(trait, "ElementValue", "Name")
                                if not name:
                                    continue
                                xref = next(_children(trait, "XRef"), None)
                                source_identifier = (
                                    xref.attrib.get("ID") if xref is not None else None
                                )
                                source_type = xref.attrib.get("DB") if xref is not None else None
                                cid = _stable(
                                    "condition", [snapshot_id, source_type, source_identifier, name]
                                )
                                conditions[cid] = ExternalCondition(
                                    condition_id=cid,
                                    snapshot_id=snapshot_id,
                                    source_identifier=source_identifier,
                                    source_identifier_type=source_type,
                                    source_name=name,
                                )
                                condition_links.append(
                                    AssertionConditionRelationship(
                                        relationship_id=_stable(
                                            "acrel", [item.assertion_instance_id, cid]
                                        ),
                                        assertion_instance_id=item.assertion_instance_id,
                                        condition_id=cid,
                                    )
                                )
                        assertions.extend(made)
                    archive.clear()
                    root.clear()
        except (ET.ParseError, EOFError, OSError) as error:
            raise ValueError("malformed ClinVar XML input") from error
        snapshot = ExternalSourceSnapshot(
            snapshot_id=snapshot_id,
            source_namespace="NCBI ClinVar",
            dataset="VCV_XML",
            release_identity=release_identity,
            release_date=release_date,
            retrieved_at=retrieved,
            source_file_sha256=exact_hash,
            source_file_size_bytes=size,
            xml_format=_local(root),
            xml_schema_identity=schema_identity,
            release_identity_override=config.release_identity_override,
        )
        variants.sort(key=lambda x: x.representation_id)
        assertions.sort(key=lambda x: x.assertion_instance_id)
        relationships.sort(key=lambda x: x.relationship_id)
        condition_links.sort(key=lambda x: x.relationship_id)
        artifacts: dict[str, str] = {}
        table_specs = {
            "external_variant_representations.parquet": (
                [x.model_dump(mode="json") for x in variants],
                {"representation_id": pl.String},
            ),
            "external_assertions.parquet": (
                [x.model_dump(mode="json") for x in assertions],
                {"assertion_instance_id": pl.String},
            ),
            "assertion_relationships.parquet": (
                [x.model_dump(mode="json") for x in relationships],
                {"relationship_id": pl.String},
            ),
            "external_conditions.parquet": (
                [x.model_dump(mode="json") for x in conditions.values()],
                {"condition_id": pl.String},
            ),
            "assertion_conditions.parquet": (
                [x.model_dump(mode="json") for x in condition_links],
                {"relationship_id": pl.String},
            ),
        }
        for name, (rows, schema) in table_specs.items():
            _write_table(temporary / name, rows, schema)
            artifacts[name] = _hash(temporary / name)
        qc = {
            "variant_representation_count": len(variants),
            "assertion_count": len(assertions),
            "assertion_levels": dict(Counter(x.assertion_level.value for x in assertions)),
            "classification_types": dict(Counter(x.classification_type.value for x in assertions)),
            "unsupported_variant_count": sum(x.variation_type != "SimpleAllele" for x in variants),
        }
        metadata = {
            "run_id": run_id,
            "snapshot": snapshot.model_dump(mode="json"),
            "parser_version": PARSER_VERSION,
            "configuration": config.model_dump(mode="json"),
            "configuration_hash": sha256(_json(config.model_dump(mode="json"))).hexdigest(),
            "package_version": __version__,
            "git_commit": _git_commit(),
            "started_at": retrieved.isoformat(),
            "completed_at": datetime.now(UTC).isoformat(),
            "xml_root_attributes": root_attrs,
        }
        report = (
            "# External evidence ingestion\n\n"
            "> Source assertions are retained claims, not project conclusions "
            "or medical advice.\n\n"
            f"- Snapshot: `{snapshot.snapshot_id}`\n- Assertions: {len(assertions)}\n"
            f"- Variant representations: {len(variants)}\n"
        )
        for name, value in (("external_source_metadata.json", metadata), ("evidence_qc.json", qc)):
            (temporary / name).write_bytes(_json(value))
            artifacts[name] = _hash(temporary / name)
        (temporary / "evidence_report.md").write_text(report)
        artifacts["evidence_report.md"] = _hash(temporary / "evidence_report.md")
        manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "snapshot_id": snapshot.snapshot_id,
            "artifacts": artifacts,
        }
        (temporary / "manifest.json").write_bytes(_json(manifest))
        temporary.rename(output_directory)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return EvidenceIngestionResult(
        run_id=run_id,
        output_directory=output_directory,
        snapshot=snapshot,
        variants=tuple(variants),
        assertions=tuple(assertions),
        relationships=tuple(relationships),
        conditions=tuple(sorted(conditions.values(), key=lambda x: x.condition_id)),
        manifest=manifest,
    )


def _make_assertion(
    classification: ET.Element,
    record: ET.Element,
    snapshot_id: str,
    representation_id: str,
    accession: str,
    version: int,
    level: AssertionLevel,
    ctype: ClassificationType,
    submitter: ET.Element | None,
    *,
    vcv_accession: str | None = None,
) -> ExternalAssertion:
    terms = _terms(classification)
    review = _first_text(classification, "ReviewStatus")
    status_term = record.attrib.get("RecordStatus") or _first_text(record, "RecordStatus")
    citations = tuple(
        sorted(
            {
                f"{n.attrib.get('Source', n.attrib.get('DB', 'unknown'))}:{n.attrib.get('ID')}"
                for n in _children(record, "ID", "XRef")
                if n.attrib.get("ID")
            }
        )
    )
    method = _first_text(record, "MethodType", "Method")
    payload = {
        "classification_element": _local(classification),
        "criteria_provided": _first_text(classification, "SubmissionName", "Citation") is not None,
    }
    fingerprint_content = {
        "accession": accession,
        "version": version,
        "type": ctype.value,
        "terms": terms,
        "review": review,
        "status": status_term,
        "submitter": dict(submitter.attrib) if submitter is not None else None,
        "date": classification.attrib.get("DateLastEvaluated"),
        "method": method,
        "citations": citations,
        "observed_in_count": sum(1 for _ in _children(record, "ObservedIn")),
    }
    fingerprint = sha256(
        json.dumps(fingerprint_content, sort_keys=True, default=str).encode()
    ).hexdigest()
    logical = f"{accession}.{version}"
    return ExternalAssertion(
        assertion_instance_id=_stable(
            "assertion", [snapshot_id, logical, fingerprint, ctype.value]
        ),
        logical_source_key=logical,
        snapshot_id=snapshot_id,
        representation_id=representation_id,
        source_accession=accession,
        accession_version=version,
        vcv_accession=vcv_accession or accession,
        scv_accession=accession if level == AssertionLevel.SUBMITTED else None,
        assertion_level=level,
        classification_type=ctype,
        source_classification_terms=terms,
        source_review_status=review,
        source_record_status=_record_status(status_term),
        source_record_status_term=status_term,
        submitter_name=(submitter.attrib.get("Name") or submitter.attrib.get("SubmitterName"))
        if submitter is not None
        else None,
        submitter_identifier=(submitter.attrib.get("ID") or submitter.attrib.get("SubmitterID"))
        if submitter is not None
        else None,
        date_last_evaluated=_parse_date(classification.attrib.get("DateLastEvaluated")),
        assertion_method=method,
        citation_identifiers=citations,
        observed_in_count=fingerprint_content["observed_in_count"],
        source_evidence_structure_count=sum(1 for _ in _children(record, "ObservedIn", "Evidence")),
        record_fingerprint=fingerprint,
        extension_payload=payload,
    )


def _validated_manifest(directory: Path, label: str) -> dict[str, Any]:
    path = directory / "manifest.json"
    if not path.is_file():
        raise ValueError(f"incomplete {label} run")
    manifest = json.loads(path.read_text())
    if not manifest.get("run_id") or not isinstance(manifest.get("artifacts"), dict):
        raise ValueError(f"invalid {label} manifest")
    for name, expected in manifest["artifacts"].items():
        artifact = directory / name
        if not artifact.is_file() or _hash(artifact) != expected:
            raise ValueError(f"{label} artifact checksum mismatch: {name}")
    return cast(dict[str, Any], manifest)


def link_external_evidence(
    m2_directory: Path,
    evidence_directory: Path,
    output_directory: Path,
    config: AnnotationConfig | None = None,
) -> AnnotationResult:
    """Link source representations to M2 variants by exact canonical allele tuple only."""
    config = config or AnnotationConfig()
    m2_directory, evidence_directory, output_directory = (
        m2_directory.resolve(),
        evidence_directory.resolve(),
        output_directory.resolve(),
    )
    m2_manifest = _validated_manifest(m2_directory, "M2")
    evidence_manifest = _validated_manifest(evidence_directory, "evidence")
    variants = pl.read_parquet(m2_directory / "variants.parquet").to_dicts()
    external = [
        ExternalVariantRepresentation.model_validate(x)
        for x in pl.read_parquet(
            evidence_directory / "external_variant_representations.parquet"
        ).to_dicts()
    ]
    index: dict[tuple[Any, ...], list[str]] = {}
    for variant in variants:
        key = (
            variant["assembly"],
            variant["chromosome"],
            variant["position"],
            variant["reference"],
            variant["alternate"],
        )
        index.setdefault(key, []).append(variant["variant_id"])
    run_id, started = str(uuid4()), datetime.now(UTC)
    temporary = _prepare_output(output_directory, run_id)
    links: list[VariantEvidenceLink] = []
    try:
        for item in sorted(external, key=lambda x: x.representation_id):
            key = (item.assembly, item.chromosome, item.position, item.reference, item.alternate)
            targets = index.get(key, [])
            if item.variation_type != "SimpleAllele":
                outcome, reason, target = (
                    LinkOutcome.UNSUPPORTED,
                    LinkReason.UNSUPPORTED_VARIATION_TYPE,
                    None,
                )
            elif None in key:
                outcome, reason, target = (
                    LinkOutcome.INCOMPATIBLE,
                    LinkReason.INCOMPLETE_ALLELE_IDENTITY,
                    None,
                )
            elif item.assembly != "GRCh38":
                outcome, reason, target = (
                    LinkOutcome.INCOMPATIBLE,
                    LinkReason.ASSEMBLY_INCOMPATIBLE,
                    None,
                )
            elif len(targets) > 1:
                outcome, reason, target = (
                    LinkOutcome.AMBIGUOUS,
                    LinkReason.MULTIPLE_CANONICAL_TARGETS,
                    None,
                )
            elif not targets:
                outcome, reason, target = (
                    LinkOutcome.UNMATCHED,
                    LinkReason.ALLELE_NOT_IN_NORMALIZATION_RUN,
                    None,
                )
            else:
                outcome, reason, target = (
                    LinkOutcome.MATCHED,
                    LinkReason.EXACT_ALLELE_MATCH,
                    targets[0],
                )
            links.append(
                VariantEvidenceLink(
                    link_id=_stable(
                        "link",
                        [
                            m2_manifest["run_id"],
                            evidence_manifest["run_id"],
                            item.representation_id,
                            key,
                            outcome.value,
                        ],
                    ),
                    annotation_run_id=run_id,
                    representation_id=item.representation_id,
                    variant_id=target,
                    outcome=outcome,
                    reason=reason,
                    assembly=item.assembly,
                    chromosome=item.chromosome,
                    position=item.position,
                    reference=item.reference,
                    alternate=item.alternate,
                )
            )
        _write_table(
            temporary / "variant_evidence_links.parquet",
            [x.model_dump(mode="json") for x in links],
            {"link_id": pl.String},
        )
        artifacts = {
            "variant_evidence_links.parquet": _hash(temporary / "variant_evidence_links.parquet")
        }
        qc = {
            "outcomes": dict(Counter(x.outcome.value for x in links)),
            "reasons": dict(Counter(x.reason.value for x in links)),
        }
        metadata = {
            "run_id": run_id,
            "m2_run_id": m2_manifest["run_id"],
            "evidence_run_id": evidence_manifest["run_id"],
            "m2_manifest_sha256": _hash(m2_directory / "manifest.json"),
            "evidence_manifest_sha256": _hash(evidence_directory / "manifest.json"),
            "algorithm_version": LINKER_VERSION,
            "configuration": config.model_dump(mode="json"),
            "configuration_hash": sha256(_json(config.model_dump(mode="json"))).hexdigest(),
            "package_version": __version__,
            "git_commit": _git_commit(),
            "source_snapshot_id": evidence_manifest.get("snapshot_id"),
            "started_at": started.isoformat(),
            "completed_at": datetime.now(UTC).isoformat(),
        }
        report = (
            "# External evidence annotation\n\n"
            "> A variant annotation does not establish that a sample carries the alternate "
            "allele, and is not a clinical conclusion.\n\n"
            + "".join(f"- {k}: {v}\n" for k, v in sorted(qc["outcomes"].items()))
        )
        for name, value in (("annotation_metadata.json", metadata), ("annotation_qc.json", qc)):
            (temporary / name).write_bytes(_json(value))
            artifacts[name] = _hash(temporary / name)
        (temporary / "annotation_report.md").write_text(report)
        artifacts["annotation_report.md"] = _hash(temporary / "annotation_report.md")
        manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "m2_run_id": m2_manifest["run_id"],
            "evidence_run_id": evidence_manifest["run_id"],
            "artifacts": artifacts,
        }
        (temporary / "manifest.json").write_bytes(_json(manifest))
        temporary.rename(output_directory)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return AnnotationResult(
        run_id=run_id, output_directory=output_directory, links=tuple(links), manifest=manifest
    )
