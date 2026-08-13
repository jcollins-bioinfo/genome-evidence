"""Offline M5 validation, exact alignment, projection, comparison, and artifacts."""

import json
import shutil
import subprocess
from collections import Counter
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
import polars as pl

from genome_evidence import __version__
from genome_evidence.domain.variants import Variant
from genome_evidence.ingest.twenty_three_and_me import _validate_private_output
from genome_evidence.normalization.models import (
    CanonicalGenotype,
    MappingCandidate,
    ObservationMapping,
)

from .models import (
    MarkerAlignment,
    MarkerAlignmentOutcome,
    MarkerExclusionReason,
    PopulationStructureConfig,
    PopulationStructureResult,
    ProjectionCoordinate,
    ProjectionDiagnostic,
    ProjectionSensitivityReplicate,
    ProjectionStatus,
    ReferenceGroupDistance,
    ReferenceNeighbor,
    ReferenceSupportEvaluation,
    ReferenceSupportStatus,
    ReferenceValidationResult,
)
from .reference import _hash, validate_population_reference

ALGORITHM_VERSION = "m5-partial-marker-lstsq-1"


def _stable(prefix: str, value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return f"{prefix}-{sha256(raw).hexdigest()}"


def _dump(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, default=str) + "\n").encode()


def _read_m2(
    directory: Path, assembly: str
) -> tuple[
    dict[str, Any], dict[str, Any], tuple[tuple[str, Variant], ...], tuple[CanonicalGenotype, ...]
]:
    manifest_path = directory / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("normalization manifest is missing")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("artifacts"), dict):
        raise ValueError("unsupported normalization schema")
    for name, expected in manifest["artifacts"].items():
        path = directory / name
        if not path.is_file() or _hash(path) != expected:
            raise ValueError(f"normalization artifact integrity failure: {name}")
    required = {
        "variants.parquet",
        "canonical_genotypes.parquet",
        "observation_mappings.parquet",
        "mapping_candidates.parquet",
        "normalization_metadata.json",
    }
    if not required <= set(manifest["artifacts"]):
        raise ValueError("normalization artifact declaration is incomplete")
    metadata = json.loads((directory / "normalization_metadata.json").read_text())
    run_id = manifest.get("run_id")
    if metadata.get("run_id") != run_id or metadata.get("target_assembly") != assembly:
        raise ValueError("normalization identity or assembly mismatch")
    try:
        variant_rows = pl.read_parquet(directory / "variants.parquet").to_dicts()
        variants = tuple(
            (str(x.pop("variant_id")), Variant.model_validate(x)) for x in variant_rows
        )
        genotypes = tuple(
            CanonicalGenotype.model_validate(x)
            for x in pl.read_parquet(directory / "canonical_genotypes.parquet").to_dicts()
        )
        mappings = tuple(
            ObservationMapping.model_validate(x)
            for x in pl.read_parquet(directory / "observation_mappings.parquet").to_dicts()
        )
        candidates = tuple(
            MappingCandidate.model_validate(x)
            for x in pl.read_parquet(directory / "mapping_candidates.parquet").to_dicts()
        )
    except Exception as error:
        raise ValueError("incompatible normalization artifact schema") from error
    ids = {v[0] for v in variants}
    if len(ids) != len(variants) or len(
        {(v.assembly, v.chromosome, v.position, v.reference, v.alternate) for _, v in variants}
    ) != len(variants):
        raise ValueError("duplicate normalization variant identity")
    all_ids = (
        [g.genotype_id for g in genotypes]
        + [m.mapping_id for m in mappings]
        + [c.candidate_id for c in candidates]
    )
    if len(all_ids) != len(set(all_ids)) or any(g.variant_id not in ids for g in genotypes):
        raise ValueError("invalid normalization scientific or genotype identity")
    mapping_ids = {m.observation_reference for m in mappings}
    candidate_ids = {c.candidate_id for c in candidates}
    if any(c.observation_reference not in mapping_ids for c in candidates) or any(
        set(m.candidate_ids) - candidate_ids for m in mappings
    ):
        raise ValueError("invalid normalization mapping relationship")
    if any(x.normalization_run_id != run_id for x in (*genotypes, *mappings, *candidates)):
        raise ValueError("cross-run normalization rows")
    return manifest, metadata, variants, genotypes


def _align(
    reference: ReferenceValidationResult,
    variants: tuple[tuple[str, Variant], ...],
    genotypes: tuple[CanonicalGenotype, ...],
) -> tuple[MarkerAlignment, ...]:
    by_key = {
        (v.assembly, v.chromosome, v.position, v.reference, v.alternate): (vid, v)
        for vid, v in variants
    }
    by_variant: dict[str, list[CanonicalGenotype]] = {}
    for genotype in genotypes:
        by_variant.setdefault(genotype.variant_id, []).append(genotype)
    aligned = []
    for marker in reference.variants:
        key = (
            marker.assembly,
            marker.chromosome,
            marker.position,
            marker.reference,
            marker.alternate,
        )
        hit = by_key.get(key)
        outcome, reason, dosage = (
            MarkerAlignmentOutcome.MISSING,
            MarkerExclusionReason.MODEL_MARKER_ABSENT,
            None,
        )
        rows: list[CanonicalGenotype] = []
        vid = None
        if hit:
            vid, variant = hit
            rows = by_variant.get(vid, [])
            valid = [
                g
                for g in rows
                if g.ploidy == 2
                and len(g.alleles) == 2
                and set(g.alleles) <= {variant.reference, variant.alternate}
            ]
            if rows and not valid:
                outcome, reason = (
                    MarkerAlignmentOutcome.EXCLUDED,
                    MarkerExclusionReason.INVALID_ALLELES,
                )
            elif valid:
                pairs = {tuple(sorted(g.alleles)) for g in valid}
                if len(pairs) > 1:
                    outcome, reason = (
                        MarkerAlignmentOutcome.DISCORDANT,
                        MarkerExclusionReason.DISCORDANT_CALLS,
                    )
                else:
                    outcome, reason = MarkerAlignmentOutcome.USED, None
                    dosage = sum(a == variant.alternate for a in valid[0].alleles)
        ids = tuple(sorted(g.genotype_id for g in rows))
        refs = tuple(sorted({g.observation_reference for g in rows}))
        aligned.append(
            MarkerAlignment(
                alignment_id=_stable(
                    "m5align",
                    [reference.identity.model_id, marker.model_variant_id, vid, outcome, ids],
                ),
                model_variant_id=marker.model_variant_id,
                variant_id=vid,
                outcome=outcome,
                exclusion_reason=reason,
                alt_dosage=dosage,
                genotype_ids=ids,
                observation_references=refs,
            )
        )
    return tuple(aligned)


def _effective(
    config: PopulationStructureConfig, reference: ReferenceValidationResult
) -> dict[str, Any]:
    m = reference.manifest
    values = {
        "components": config.component_count or m.default_component_count,
        "marker_count": max(
            config.minimum_observed_marker_count or 0, m.minimum_observed_marker_count
        ),
        "marker_fraction": max(
            config.minimum_observed_marker_fraction or 0, m.minimum_observed_marker_fraction
        ),
        "chromosomes": max(config.minimum_chromosomes or 0, m.minimum_chromosomes),
        "loading_energy": max(config.minimum_loading_energy or 0, m.minimum_loading_energy),
        "condition": min(
            config.maximum_condition_number or m.maximum_condition_number,
            m.maximum_condition_number,
        ),
        "sensitivity": max(
            config.sensitivity_minimum_valid_replicates or 0, m.sensitivity_minimum_valid_replicates
        ),
    }
    if values["components"] > m.maximum_component_count:
        raise ValueError("requested component count exceeds reference limit")
    return values


def _project(
    reference: ReferenceValidationResult,
    alignments: tuple[MarkerAlignment, ...],
    thresholds: dict[str, Any],
    omit: str | None = None,
) -> tuple[ProjectionStatus, ProjectionDiagnostic, np.ndarray | None]:
    k = thresholds["components"]
    by_id = {x.model_variant_id: x for x in alignments}
    used = [
        (v, by_id[v.model_variant_id])
        for v in reference.variants
        if by_id[v.model_variant_id].outcome == MarkerAlignmentOutcome.USED and v.chromosome != omit
    ]
    matrix = np.asarray([v.loadings[:k] for v, _ in used], dtype=np.float64)
    x = np.asarray(
        [(a.alt_dosage - v.training_mean_dosage) / v.training_scale for v, a in used],
        dtype=np.float64,
    )
    all_l = np.asarray([v.loadings[:k] for v in reference.variants], dtype=np.float64)
    energy = (
        tuple((np.square(matrix).sum(axis=0) / np.square(all_l).sum(axis=0)).tolist())
        if used
        else tuple(0.0 for _ in range(k))
    )
    chromosomes = tuple(sorted({v.chromosome for v, _ in used}, key=int))
    rank = int(np.linalg.matrix_rank(matrix)) if used else 0
    singular = tuple(float(x) for x in np.linalg.svd(matrix, compute_uv=False)) if used else ()
    condition = float(np.linalg.cond(matrix)) if used and rank == k else None
    failed = []
    fraction = len(used) / len(reference.variants)
    if len(used) < thresholds["marker_count"]:
        failed.append("minimum_observed_marker_count")
    if fraction < thresholds["marker_fraction"]:
        failed.append("minimum_observed_marker_fraction")
    if len(chromosomes) < thresholds["chromosomes"]:
        failed.append("minimum_chromosomes")
    if any(e < thresholds["loading_energy"] for e in energy):
        failed.append("minimum_loading_energy")
    if rank < k:
        failed.append("full_rank")
    if condition is not None and condition > thresholds["condition"]:
        failed.append("maximum_condition_number")
    if not np.isfinite(matrix).all() or not np.isfinite(x).all():
        failed.append("finite_inputs")
    coordinates = None
    residual = None
    if not failed:
        coordinates, residuals, _, _ = np.linalg.lstsq(
            matrix, x, rcond=reference.manifest.least_squares_rcond
        )
        residual = float(np.linalg.norm(matrix @ coordinates - x))
        if not np.isfinite(coordinates).all() or not np.isfinite(residual):
            failed.append("finite_outputs")
            coordinates = None
    diagnostic = ProjectionDiagnostic(
        selected_component_count=k,
        used_marker_count=len(used),
        observed_marker_fraction=fraction,
        chromosomes=chromosomes,
        loading_energy_coverage=energy,
        matrix_rank=rank,
        singular_values=singular,
        condition_number=condition,
        residual_norm=residual,
        failed_gates=tuple(failed),
    )
    return (
        (ProjectionStatus.NOT_PROJECTED if failed else ProjectionStatus.PROJECTED),
        diagnostic,
        coordinates,
    )


def _compare(
    reference: ReferenceValidationResult, z: np.ndarray, count: int
) -> tuple[
    tuple[ReferenceGroupDistance, ...],
    tuple[ReferenceNeighbor, ...],
    tuple[ReferenceSupportEvaluation, ...],
    ReferenceSupportStatus,
]:
    k = len(z)
    scores = np.asarray([x[2][:k] for x in reference.reference_scores])
    sd = np.asarray([x.reference_score_standard_deviation for x in reference.components[:k]])
    groups = {g.reference_group_id: g for g in reference.groups}
    centroids = {
        gid: scores[[x[1] == gid for x in reference.reference_scores]].mean(axis=0)
        for gid in groups
    }
    raw_group = sorted(
        ((float(np.linalg.norm((z - c) / sd)), gid) for gid, c in centroids.items()),
        key=lambda x: (x[0], x[1]),
    )
    distances = tuple(
        ReferenceGroupDistance(
            reference_group_id=gid,
            source_label=groups[gid].source_label,
            source_definition=groups[gid].source_definition,
            label_provenance=groups[gid].label_provenance,
            sample_count=groups[gid].declared_sample_count,
            distance=d,
            rank=i,
        )
        for i, (d, gid) in enumerate(raw_group, 1)
    )
    raw_samples = sorted(
        (
            (float(np.linalg.norm((z - s) / sd)), sid, gid)
            for s, (sid, gid, _) in zip(scores, reference.reference_scores, strict=True)
        ),
        key=lambda x: (x[0], x[2], x[1]),
    )[:count]
    neighbors = tuple(
        ReferenceNeighbor(reference_sample_id=sid, reference_group_id=gid, distance=d, rank=i)
        for i, (d, sid, gid) in enumerate(raw_samples, 1)
    )
    support = []
    for gid, _group in sorted(groups.items()):
        members = scores[[x[1] == gid for x in reference.reference_scores]]
        d = float(np.linalg.norm((z - centroids[gid]) / sd))
        if len(members) < reference.manifest.support_minimum_group_size:
            support.append(
                ReferenceSupportEvaluation(
                    reference_group_id=gid,
                    evaluable=False,
                    inside_envelope=None,
                    distance=None,
                    envelope_radius=None,
                )
            )
            continue
        radii = np.linalg.norm((members - centroids[gid]) / sd, axis=1)
        radius = float(
            np.quantile(
                radii,
                reference.manifest.support_quantile,
                method=reference.manifest.support_quantile_method,
            )
        )
        support.append(
            ReferenceSupportEvaluation(
                reference_group_id=gid,
                evaluable=True,
                inside_envelope=d <= radius,
                distance=d,
                envelope_radius=radius,
            )
        )
    evaluable = [s for s in support if s.evaluable]
    status = (
        ReferenceSupportStatus.NOT_EVALUABLE
        if not evaluable
        else ReferenceSupportStatus.WITHIN
        if any(s.inside_envelope for s in evaluable)
        else ReferenceSupportStatus.OUTSIDE
    )
    return distances, neighbors, tuple(support), status


def infer_population_structure(
    normalization_directory: Path,
    reference_bundle_directory: Path,
    output_directory: Path,
    config: PopulationStructureConfig | None = None,
) -> PopulationStructureResult:
    config = config or PopulationStructureConfig()
    _validate_private_output(output_directory)
    if output_directory.exists():
        raise FileExistsError("population-structure output already exists")
    reference = validate_population_reference(reference_bundle_directory)
    m2, metadata, variants, genotypes = _read_m2(
        normalization_directory, reference.manifest.assembly
    )
    thresholds = _effective(config, reference)
    alignments = _align(reference, variants, genotypes)
    status, diagnostic, z = _project(reference, alignments, thresholds)
    coordinates = ()
    distances = ()
    neighbors = ()
    support = ()
    support_status = None
    sensitivity = []
    if z is not None:
        coordinates = tuple(
            ProjectionCoordinate(
                coordinate_id=_stable(
                    "m5pc",
                    [
                        reference.identity.model_id,
                        i,
                        float(value),
                        [
                            a.alignment_id
                            for a in alignments
                            if a.outcome == MarkerAlignmentOutcome.USED
                        ],
                    ],
                ),
                component_index=i,
                coordinate=float(value),
            )
            for i, value in enumerate(z, 1)
        )
        distances, neighbors, support, support_status = _compare(
            reference, z, config.nearest_neighbor_count
        )
        for chromosome in diagnostic.chromosomes:
            rs, rd, rz = _project(reference, alignments, thresholds, chromosome)
            rg = None
            rss = None
            if rz is not None:
                gd, _, _, rss = _compare(reference, rz, config.nearest_neighbor_count)
                rg = gd[0].reference_group_id
            sensitivity.append(
                ProjectionSensitivityReplicate(
                    omitted_chromosome=chromosome,
                    status=rs,
                    retained_marker_count=rd.used_marker_count,
                    retained_marker_fraction=rd.observed_marker_fraction,
                    diagnostic=rd,
                    coordinates=tuple(float(v) for v in rz) if rz is not None else (),
                    absolute_deviations=tuple(float(v) for v in np.abs(rz - z))
                    if rz is not None
                    else (),
                    rank_one_group_id=rg,
                    support_status=rss,
                )
            )
    run = str(uuid4())
    started = datetime.now(UTC)
    temp = output_directory.with_name(f".{output_directory.name}.{run}.tmp")
    temp.mkdir(parents=True)
    try:
        tables = {
            "marker_alignment.parquet": [x.model_dump(mode="json") for x in alignments],
            "projected_coordinates.parquet": [x.model_dump(mode="json") for x in coordinates],
            "reference_group_distances.parquet": [x.model_dump(mode="json") for x in distances],
            "reference_neighbors.parquet": [x.model_dump(mode="json") for x in neighbors],
            "reference_support.parquet": [x.model_dump(mode="json") for x in support],
            "projection_sensitivity.parquet": [x.model_dump(mode="json") for x in sensitivity],
        }
        schemas = {name: pl.Schema({"empty": pl.String}) for name in tables}
        hashes = {}
        for name, rows in tables.items():
            (pl.DataFrame(rows) if rows else pl.DataFrame(schema=schemas[name])).write_parquet(
                temp / name
            )
            hashes[name] = {"sha256": _hash(temp / name), "byte_size": (temp / name).stat().st_size}
        qc = {
            "projection_status": status,
            "support_status": support_status,
            "alignment_outcomes": dict(Counter(x.outcome.value for x in alignments)),
            "diagnostic": diagnostic.model_dump(mode="json"),
            "sensitivity_eligible": len(sensitivity),
            "sensitivity_valid": sum(x.status == ProjectionStatus.PROJECTED for x in sensitivity),
            "sensitivity_summary_available": sum(
                x.status == ProjectionStatus.PROJECTED for x in sensitivity
            )
            >= thresholds["sensitivity"],
        }
        warning = (
            "Results depend on the exact supplied panel, marker selection, preprocessing, "
            "labels, PCA axes, and model version. Reference labels describe sampled subsets, "
            "not discrete natural kinds. Genetic similarity is not race, ethnicity, nationality, "
            "culture, or identity. A nearest reference subset is not an ancestry assignment. "
            "No ancestry percentages were estimated. Missing array markers and SNP ascertainment "
            "can distort projection. PCA can reflect technical artifacts, relatedness, sample "
            "imbalance, LD, and outliers. Projected PCs can shrink or become unstable, especially "
            "for later components. Outside-support means only poor representation by this model. "
            "This is not a diagnosis, disease-risk estimate, medical recommendation, or negative "
            "genetic test."
        )
        report = "# Reference-panel population-structure projection\n\n> " + warning + "\n\n"
        report += f"- Projection status: `{status}`\n"
        report += f"- Support status: `{support_status}`\n"
        report += f"- Used markers: {diagnostic.used_marker_count}\n"
        report += f"- Failed gates: {', '.join(diagnostic.failed_gates) or 'none'}\n"
        config_hash = sha256(
            json.dumps(
                config.model_dump(mode="json"), sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest()
        md = {
            "run_id": run,
            "m2_run_id": m2["run_id"],
            "m2_manifest_sha256": _hash(normalization_directory / "manifest.json"),
            "m2_artifacts": m2["artifacts"],
            "reference_identity": reference.identity.model_dump(mode="json"),
            "configuration": config.model_dump(mode="json"),
            "configuration_hash": config_hash,
            "algorithm_version": ALGORITHM_VERSION,
            "package_version": __version__,
            "git_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], capture_output=True, text=True
            ).stdout.strip(),
            "started_at": started.isoformat(),
            "completed_at": datetime.now(UTC).isoformat(),
            "projection_status": status,
            "support_status": support_status,
            "diagnostic": diagnostic.model_dump(mode="json"),
        }
        for name, data in {
            "population_structure_qc.json": _dump(qc),
            "population_structure_metadata.json": _dump(md),
            "population_structure_report.md": report.encode(),
        }.items():
            (temp / name).write_bytes(data)
            hashes[name] = {"sha256": _hash(temp / name), "byte_size": (temp / name).stat().st_size}
        manifest = {
            "schema_version": 1,
            "run_id": run,
            "m2_run_id": m2["run_id"],
            "reference_model_id": reference.identity.model_id,
            "reference_model_version": reference.identity.model_version,
            "configuration_hash": config_hash,
            "artifacts": hashes,
        }
        (temp / "manifest.json").write_bytes(_dump(manifest))
        temp.rename(output_directory)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    return PopulationStructureResult(
        run_id=run,
        output_directory=output_directory,
        reference_identity=reference.identity,
        projection_status=status,
        support_status=support_status,
        alignments=alignments,
        diagnostic=diagnostic,
        coordinates=coordinates,
        group_distances=distances,
        neighbors=neighbors,
        support=support,
        sensitivity=tuple(sensitivity),
    )
