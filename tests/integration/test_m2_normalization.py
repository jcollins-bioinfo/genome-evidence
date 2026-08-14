import json
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from genome_evidence.cli import app
from genome_evidence.ingest import Ingest23andMeConfig, ingest_23andme
from genome_evidence.normalization import NormalizationConfig, normalize_m1_run


def resources(path: Path) -> tuple[Path, Path]:
    markers = path / "synthetic-markers.json"
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
                {
                    "marker_id": "internalSynthetic",
                    "assembly": "GRCh38",
                    "chromosome": "X",
                    "position": 303,
                    "reference": "A",
                    "alternate": "C",
                    "orientation": "none",
                    "orientation_authoritative": True,
                },
            ]
        )
    )
    fasta = path / "synthetic.fa"
    fasta.write_text(
        ">1\n" + "A" * 100 + "A" + "A" * 100 + "C" + "A" * 200 + "\n>X\n" + "A" * 400 + "\n"
    )
    return markers, fasta


def test_end_to_end_same_build_and_no_call(tmp_path: Path) -> None:
    fixture = Path(__file__).parents[1] / "fixtures" / "23andme" / "normal.txt"
    m1 = tmp_path / "m1"
    ingest_23andme(fixture, m1, Ingest23andMeConfig(genome_build_override="GRCh38"))
    markers, fasta = resources(tmp_path)
    result = normalize_m1_run(
        m1, tmp_path / "m2", NormalizationConfig(marker_definitions=markers, target_reference=fasta)
    )
    assert len(result.mappings) == 3
    assert len({x.observation_reference for x in result.mappings}) == 3
    assert len(result.genotypes) == 2
    assert all(
        "raw_genotype" not in (tmp_path / "m2" / "normalization_report.md").read_text() for _ in [0]
    )
    assert pl.read_parquet(tmp_path / "m2" / "observation_mappings.parquet").height == 3
    manifest = json.loads((tmp_path / "m2" / "manifest.json").read_text())
    from hashlib import sha256

    for name, digest in manifest["artifacts"].items():
        assert sha256((tmp_path / "m2" / name).read_bytes()).hexdigest() == digest


def test_checksum_rejected_and_cli(tmp_path: Path) -> None:
    fixture = Path(__file__).parents[1] / "fixtures" / "23andme" / "normal.txt"
    m1 = tmp_path / "m1"
    ingest_23andme(fixture, m1)
    markers, fasta = resources(tmp_path)
    (m1 / "observations.parquet").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="checksum"):
        normalize_m1_run(
            m1,
            tmp_path / "m2",
            NormalizationConfig(marker_definitions=markers, target_reference=fasta),
        )
    cli = CliRunner().invoke(
        app,
        [
            "normalize",
            "--input",
            str(m1),
            "--output",
            str(tmp_path / "cli"),
            "--marker-definitions",
            str(markers),
            "--target-reference",
            str(fasta),
        ],
    )
    assert cli.exit_code != 0


def test_duplicate_observation_references_use_source_line(tmp_path: Path) -> None:
    fixture = Path(__file__).parents[1] / "fixtures" / "23andme" / "duplicate_marker.txt"
    m1 = tmp_path / "m1"
    ingest_23andme(fixture, m1, Ingest23andMeConfig(genome_build_override="GRCh38"))
    markers = tmp_path / "m.json"
    markers.write_text(
        json.dumps(
            [
                {
                    "marker_id": "syntheticDup",
                    "assembly": "GRCh38",
                    "chromosome": "1",
                    "position": 10,
                    "reference": "A",
                    "alternate": "G",
                    "orientation": "none",
                    "orientation_authoritative": True,
                }
            ]
        )
    )
    fasta = tmp_path / "r.fa"
    fasta.write_text(">1\n" + "A" * 20 + "\n")
    first = normalize_m1_run(
        m1,
        tmp_path / "one",
        NormalizationConfig(marker_definitions=markers, target_reference=fasta),
    )
    second = normalize_m1_run(
        m1,
        tmp_path / "two",
        NormalizationConfig(marker_definitions=markers, target_reference=fasta),
    )
    assert first.mappings[0].observation_reference != first.mappings[1].observation_reference
    assert [x.observation_reference for x in first.mappings] == [
        x.observation_reference for x in second.mappings
    ]


def test_marker_definition_must_match_resolved_source_assembly(tmp_path: Path) -> None:
    fixture = Path(__file__).parents[1] / "fixtures" / "23andme" / "duplicate_marker.txt"
    m1 = tmp_path / "m1"
    ingest_23andme(fixture, m1, Ingest23andMeConfig(genome_build_override="GRCh38"))
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            [
                {
                    "marker_id": "syntheticDup",
                    "assembly": "GRCh37",
                    "chromosome": "1",
                    "position": 10,
                    "reference": "A",
                    "alternate": "G",
                    "orientation": "none",
                    "orientation_authoritative": True,
                }
            ]
        )
    )
    fasta = tmp_path / "GRCh38.fa"
    fasta.write_text(">1\n" + "A" * 20 + "\n")

    result = normalize_m1_run(
        m1,
        tmp_path / "m2",
        NormalizationConfig(marker_definitions=markers, target_reference=fasta),
    )

    assert {mapping.outcome.value for mapping in result.mappings} == {"unsupported"}
    assert {mapping.reason for mapping in result.mappings} == {
        "MARKER_DEFINITION_ASSEMBLY_MISMATCH"
    }
