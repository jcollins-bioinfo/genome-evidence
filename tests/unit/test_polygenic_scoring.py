import json
from hashlib import sha256
from pathlib import Path

import polars as pl

from genome_evidence.ingest import Ingest23andMeConfig, ingest_23andme
from genome_evidence.normalization import NormalizationConfig, normalize_m1_run
from genome_evidence.polygenic_scoring import (
    ScoreConfig,
    calculate_polygenic_scores,
    validate_polygenic_score_bundle,
)
from genome_evidence.polygenic_scoring.reference import compare_reference


def _bundle(root: Path) -> Path:
    root.mkdir()
    metadata = {
        "models": [
            {
                "pgs_id": "PGS999999",
                "version": "1",
                "trait": "fabricated trait",
                "source_url": "https://example.invalid/fabricated",
                "citation": "synthetic fixture; no publication",
                "license": "CC0 synthetic fixture",
                "assembly": "GRCh38",
                "declared_variant_count": 3,
                "completeness_policy": "all_supported_variants_required",
            }
        ]
    }
    (root / "model_metadata.json").write_text(json.dumps(metadata))
    pl.DataFrame(
        [
            {
                "pgs_id": "PGS999999",
                "assembly": "GRCh38",
                "chromosome": "1",
                "position": 101,
                "reference": "A",
                "alternate": "G",
                "effect_allele": "G",
                "effect_weight": "1.5",
            },
            {
                "pgs_id": "PGS999999",
                "assembly": "GRCh38",
                "chromosome": "1",
                "position": 202,
                "reference": "C",
                "alternate": "T",
                "effect_allele": "C",
                "effect_weight": "-2.5e-1",
            },
            {
                "pgs_id": "PGS999999",
                "assembly": "GRCh38",
                "chromosome": "2",
                "position": 10,
                "reference": "A",
                "alternate": "C",
                "effect_allele": "C",
                "effect_weight": "0.25",
            },
        ]
    ).write_parquet(root / "model_variants.parquet")
    pl.DataFrame(schema={"pgs_id": pl.String, "reason": pl.String}).write_parquet(
        root / "model_exclusions.parquet"
    )
    artifacts = {}
    for name in ("model_metadata.json", "model_variants.parquet", "model_exclusions.parquet"):
        path = root / name
        artifacts[name] = {
            "byte_size": path.stat().st_size,
            "sha256": sha256(path.read_bytes()).hexdigest(),
            "privacy_class": "public_or_controlled_aggregate",
            "redistribution": "synthetic_fixture",
        }
    (root / "manifest.json").write_text(
        json.dumps(
            {
                "schema": "genome-evidence-pgs-bundle/v1",
                "bundle_id": "synthetic-pgs-bundle-v1",
                "artifacts": artifacts,
            }
        )
    )
    return root


def _m2(root: Path) -> Path:
    fixture = Path("tests/fixtures/23andme/normal.txt")
    m1 = root / "m1"
    ingest_23andme(fixture, m1, Ingest23andMeConfig(genome_build_override="GRCh38"))
    markers = root / "markers.json"
    markers.write_text(
        json.dumps(
            [
                {
                    "marker_id": "rsSynthetic1",
                    "assembly": "GRCh38",
                    "chromosome": "1",
                    "position": 101,
                    "reference": "A",
                    "alternate": "G",
                    "orientation": "none",
                    "orientation_authoritative": True,
                },
                {
                    "marker_id": "rsSynthetic2",
                    "assembly": "GRCh38",
                    "chromosome": "1",
                    "position": 202,
                    "reference": "C",
                    "alternate": "T",
                    "orientation": "none",
                    "orientation_authoritative": True,
                },
            ]
        )
    )
    fasta = root / "ref.fa"
    fasta.write_text(">1\n" + "A" * 100 + "A" + "A" * 100 + "C" + "A" * 20 + "\n")
    normalize_m1_run(
        m1, root / "m2", NormalizationConfig(marker_definitions=markers, target_reference=fasta)
    )
    return root / "m2"


def test_exact_scoring_ref_orientation_missing_and_completion(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path / "bundle")
    m2 = _m2(tmp_path)
    assert validate_polygenic_score_bundle(bundle).model_ids == ("PGS999999",)
    calculate_polygenic_scores(m2, bundle, tmp_path / "out", ScoreConfig(pgs_ids=("PGS999999",)))
    row = pl.read_parquet(tmp_path / "out/polygenic_score_results.parquet").to_dicts()[0]
    assert row["raw_partial_score"] == "1"
    assert row["status"] == "partial_not_evaluable"
    assert row["contributed_marker_count"] == 2
    assert pl.read_parquet(tmp_path / "out/score_exclusions.parquet")["reason"].to_list() == [
        "TARGET_VARIANT_MISSING"
    ]
    assert (tmp_path / "out/COMPLETED.json").is_file()


def test_reference_requires_exact_compatibility() -> None:
    identity = {
        "pgs_id": "PGS999999",
        "model_version": "1",
        "assembly": "GRCh38",
        "matching_pipeline": "exact-v1",
        "missingness_policy": "complete",
    }
    assert (
        compare_reference(2, {**identity, "scores": [0, 1, 2, 3]}, identity)["status"]
        == "evaluable"
    )
    assert (
        compare_reference(2, {**identity, "assembly": "GRCh37", "scores": [0, 1]}, identity)[
            "status"
        ]
        == "not_evaluable"
    )
    assert (
        compare_reference(2, {**identity, "scores": [1, 1]}, identity)["reason"]
        == "REFERENCE_ZERO_VARIANCE"
    )
