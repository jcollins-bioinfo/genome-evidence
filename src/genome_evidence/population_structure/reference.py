"""Strict validation for offline, checksummed PCA reference bundles."""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import polars as pl
from pydantic import ValidationError

from .models import (
    ComponentMetadata,
    PopulationReferenceIdentity,
    ReferenceBundleManifest,
    ReferenceGroup,
    ReferenceValidationResult,
    ReferenceVariantLoading,
)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_population_reference(reference_bundle_directory: Path) -> ReferenceValidationResult:
    manifest_path = reference_bundle_directory / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("reference manifest is missing")
    try:
        manifest = ReferenceBundleManifest.model_validate_json(manifest_path.read_bytes())
    except (ValidationError, json.JSONDecodeError) as error:
        raise ValueError("invalid reference manifest") from error
    for name, identity in manifest.artifacts.items():
        path = reference_bundle_directory / name
        if (
            not path.is_file()
            or path.stat().st_size != identity.byte_size
            or _hash(path) != identity.sha256
        ):
            raise ValueError(f"reference artifact integrity failure: {name}")
    try:
        variants = tuple(
            ReferenceVariantLoading.model_validate(x)
            for x in pl.read_parquet(
                reference_bundle_directory / "variant_loadings.parquet"
            ).to_dicts()
        )
        groups = tuple(
            ReferenceGroup.model_validate(x)
            for x in pl.read_parquet(
                reference_bundle_directory / "reference_groups.parquet"
            ).to_dicts()
        )
        components = tuple(
            ComponentMetadata.model_validate(x)
            for x in pl.read_parquet(
                reference_bundle_directory / "component_metadata.parquet"
            ).to_dicts()
        )
        score_rows = pl.read_parquet(
            reference_bundle_directory / "reference_scores.parquet"
        ).to_dicts()
        scores = tuple(
            (
                str(x["reference_sample_id"]),
                str(x["reference_group_id"]),
                tuple(float(v) for v in x["scores"]),
            )
            for x in score_rows
        )
    except Exception as error:
        raise ValueError("invalid reference artifact schema") from error
    if (len(variants), len(scores), len(groups), len(components)) != (
        manifest.marker_count,
        manifest.reference_sample_count,
        manifest.reference_group_count,
        manifest.component_count,
    ):
        raise ValueError("reference artifact counts disagree with manifest")
    keys = [(v.assembly, v.chromosome, v.position, v.reference, v.alternate) for v in variants]
    if len({v.model_variant_id for v in variants}) != len(variants) or len(set(keys)) != len(keys):
        raise ValueError("duplicate reference marker identity")
    if len({x[0] for x in scores}) != len(scores) or len(
        {g.reference_group_id for g in groups}
    ) != len(groups):
        raise ValueError("duplicate reference sample or group identity")
    group_ids = {g.reference_group_id for g in groups}
    if any(s[1] not in group_ids for s in scores):
        raise ValueError("unknown reference group membership")
    if any(
        sum(s[1] == g.reference_group_id for s in scores) != g.declared_sample_count for g in groups
    ):
        raise ValueError("reference group sample count mismatch")
    bases = set("ACGT")
    if any(
        v.assembly != "GRCh38"
        or v.chromosome not in {str(i) for i in range(1, 23)}
        or len(v.reference) != 1
        or len(v.alternate) != 1
        or v.reference not in bases
        or v.alternate not in bases
        or v.reference == v.alternate
        or v.effect_allele != v.alternate
        for v in variants
    ):
        raise ValueError("unsupported reference marker")
    k = manifest.component_count
    if (
        [c.component_index for c in components] != list(range(1, k + 1))
        or any(len(v.loadings) != k for v in variants)
        or any(len(s[2]) != k for s in scores)
    ):
        raise ValueError("reference component dimensions are inconsistent")
    arrays = [
        np.asarray([v.loadings for v in variants], dtype=np.float64),
        np.asarray([s[2] for s in scores], dtype=np.float64),
    ]
    if any(not np.isfinite(a).all() for a in arrays):
        raise ValueError("reference contains nonfinite values")
    loadings, score_matrix = arrays
    if not np.allclose(loadings.T @ loadings, np.eye(k), atol=manifest.numerical_tolerance, rtol=0):
        raise ValueError("loading columns are not orthonormal")
    means, sds = score_matrix.mean(axis=0), score_matrix.std(axis=0)
    expected_means = np.asarray([c.reference_score_mean for c in components])
    expected_sds = np.asarray([c.reference_score_standard_deviation for c in components])
    if not np.allclose(
        means, expected_means, atol=manifest.numerical_tolerance, rtol=0
    ) or not np.allclose(sds, expected_sds, atol=manifest.numerical_tolerance, rtol=0):
        raise ValueError("reference score summaries disagree with component metadata")
    identity = PopulationReferenceIdentity(
        model_id=manifest.model_id,
        model_version=manifest.model_version,
        assembly=manifest.assembly,
        manifest_sha256=_hash(manifest_path),
        artifacts=manifest.artifacts,
    )
    return ReferenceValidationResult(
        directory=reference_bundle_directory,
        identity=identity,
        manifest=manifest,
        variants=variants,
        groups=groups,
        components=components,
        reference_scores=scores,
    )
