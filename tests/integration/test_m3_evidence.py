import gzip
import json
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from genome_evidence.cli import app
from genome_evidence.evidence import ingest_clinvar_vcv, link_external_evidence
from genome_evidence.evidence.models import AssertionLevel, ClassificationType, LinkOutcome
from genome_evidence.ingest import Ingest23andMeConfig, ingest_23andme
from genome_evidence.normalization import NormalizationConfig, normalize_m1_run

FIXTURE = Path(__file__).parents[1] / "fixtures" / "clinvar" / "synthetic_vcv.xml"


def _m2(tmp_path: Path) -> Path:
    source = tmp_path / "synthetic.txt"
    source.write_text("# genome build: GRCh38\nsynthetic-evidence-marker\t1\t101\tAA\n")
    m1 = tmp_path / "m1"
    ingest_23andme(source, m1, Ingest23andMeConfig(genome_build_override="GRCh38"))
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            [
                {
                    "marker_id": "synthetic-evidence-marker",
                    "assembly": "GRCh38",
                    "chromosome": "1",
                    "position": 101,
                    "reference": "A",
                    "alternate": "G",
                    "orientation": "none",
                    "orientation_authoritative": True,
                }
            ]
        )
    )
    fasta = tmp_path / "synthetic.fa"
    fasta.write_text(">1\n" + "A" * 200 + "\n")
    m2 = tmp_path / "m2"
    normalize_m1_run(
        m1, m2, NormalizationConfig(marker_definitions=markers, target_reference=fasta)
    )
    return m2


def test_plain_gzip_assertions_and_deterministic_scientific_ids(tmp_path: Path) -> None:
    plain = ingest_clinvar_vcv(FIXTURE, tmp_path / "plain")
    gz = tmp_path / "synthetic.xml.gz"
    with gzip.GzipFile(filename=str(gz), mode="wb", mtime=0) as stream:
        stream.write(FIXTURE.read_bytes())
    compressed = ingest_clinvar_vcv(gz, tmp_path / "gz")
    assert len(plain.variants) == 2
    assert {a.assertion_level for a in plain.assertions} == {
        AssertionLevel.AGGREGATE,
        AssertionLevel.SUBMITTED,
    }
    assert {a.classification_type for a in plain.assertions} >= {
        ClassificationType.GERMLINE,
        ClassificationType.SOMATIC_CLINICAL_IMPACT,
        ClassificationType.ONCOGENICITY,
    }
    scvs = [a for a in plain.assertions if a.scv_accession]
    assert {a.logical_source_key for a in scvs} >= {"SCV999000001.1", "SCV999000002.3"}
    assert len({a.assertion_instance_id for a in scvs}) == len(scvs)
    assert (
        plain.snapshot.snapshot_id != compressed.snapshot.snapshot_id
    )  # exact supplied bytes differ
    assert any(a.submitter_name == "Synthetic Lab Alpha" for a in scvs)
    assert plain.relationships


def test_exact_link_is_not_carriage_and_checksums_are_validated(tmp_path: Path) -> None:
    evidence_dir = tmp_path / "evidence"
    ingest_clinvar_vcv(FIXTURE, evidence_dir)
    m2 = _m2(tmp_path)
    annotation = link_external_evidence(m2, evidence_dir, tmp_path / "annotation")
    assert [link.outcome for link in annotation.links] == [
        LinkOutcome.MATCHED,
        LinkOutcome.UNSUPPORTED,
    ]
    assert (
        "does not establish that a sample carries"
        in (tmp_path / "annotation" / "annotation_report.md").read_text()
    )
    assert (
        "genotype"
        not in pl.read_parquet(tmp_path / "annotation" / "variant_evidence_links.parquet").columns
    )
    (evidence_dir / "external_assertions.parquet").write_bytes(b"tampered")
    with pytest.raises(ValueError, match="checksum"):
        link_external_evidence(m2, evidence_dir, tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()


def test_wrong_root_malformed_existing_output_and_cli(tmp_path: Path) -> None:
    wrong = tmp_path / "wrong.xml"
    wrong.write_text("<NotClinVar/>")
    with pytest.raises(ValueError, match="unsupported"):
        ingest_clinvar_vcv(wrong, tmp_path / "wrong-run")
    malformed = tmp_path / "malformed.xml"
    malformed.write_text("<ReleaseSet")
    with pytest.raises(ValueError, match="malformed"):
        ingest_clinvar_vcv(malformed, tmp_path / "malformed-run")
    output = tmp_path / "cli"
    result = CliRunner().invoke(
        app, ["evidence", "ingest-clinvar", "--input", str(FIXTURE), "--output", str(output)]
    )
    assert result.exit_code == 0, result.output
    again = CliRunner().invoke(
        app, ["evidence", "ingest-clinvar", "--input", str(FIXTURE), "--output", str(output)]
    )
    assert again.exit_code == 2
