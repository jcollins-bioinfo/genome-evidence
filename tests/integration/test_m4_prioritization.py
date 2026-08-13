import json
from hashlib import sha256
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from genome_evidence.cli import app
from genome_evidence.evidence import ingest_clinvar_vcv, link_external_evidence
from genome_evidence.ingest import Ingest23andMeConfig, ingest_23andme
from genome_evidence.normalization import NormalizationConfig, normalize_m1_run
from genome_evidence.prioritization import prioritize_clinical_variants
from genome_evidence.prioritization.models import (
    AnalysisContext,
    CandidateEligibility,
    ClinicalPrioritizationConfig,
    GenotypeEvidenceState,
    PrioritizationResult,
    ReviewPriorityBand,
)

POLICY = Path(__file__).parents[2] / "references/clinvar-germline-review-policy-v1.json"


def _runs(tmp_path: Path) -> tuple[Path, Path, Path]:
    calls = [
        ("alt", 10, "AG"),
        ("ref", 20, "CC"),
        ("missing", 30, "--"),
        ("discord", 40, "AG"),
        ("discord", 40, "AA"),
        ("benign", 50, "AT"),
        ("somatic", 60, "AC"),
    ]
    source = tmp_path / "calls.txt"
    source.write_text(
        "# genome build: GRCh38\n"
        + "".join(f"{marker}\t1\t{position}\t{call}\n" for marker, position, call in calls)
    )
    markers = tmp_path / "markers.json"
    markers.write_text(
        json.dumps(
            [
                {
                    "marker_id": marker,
                    "assembly": "GRCh38",
                    "chromosome": "1",
                    "position": pos,
                    "reference": ref,
                    "alternate": alt,
                    "orientation": "none",
                    "orientation_authoritative": True,
                }
                for marker, pos, ref, alt in [
                    ("alt", 10, "A", "G"),
                    ("ref", 20, "C", "T"),
                    ("missing", 30, "G", "A"),
                    ("discord", 40, "A", "G"),
                    ("benign", 50, "A", "T"),
                    ("somatic", 60, "A", "C"),
                ]
            ]
        )
    )
    fasta = tmp_path / "ref.fa"
    sequence = list("A" * 100)
    sequence[19], sequence[29] = "C", "G"
    fasta.write_text(">1\n" + "".join(sequence) + "\n")
    ingest_23andme(source, tmp_path / "m1", Ingest23andMeConfig(genome_build_override="GRCh38"))
    normalize_m1_run(
        tmp_path / "m1",
        tmp_path / "m2",
        NormalizationConfig(marker_definitions=markers, target_reference=fasta),
    )
    records = [
        (1, 10, "A", "G", "Pathogenic", "GermlineClassification", "current"),
        (2, 20, "C", "T", "Pathogenic", "GermlineClassification", "current"),
        (3, 30, "G", "A", "Pathogenic", "GermlineClassification", "current"),
        (4, 40, "A", "G", "Uncertain significance", "GermlineClassification", "current"),
        (5, 50, "A", "T", "Benign", "GermlineClassification", "current"),
        (6, 60, "A", "C", "Tier I", "SomaticClinicalImpact", "current"),
    ]
    xml = tmp_path / "clinvar.xml"
    body = []
    for number, pos, ref, alt, term, classification, status in records:
        body.append(
            f'<VariationArchive Accession="VCV90000000{number}" Version="1" '
            f'RecordStatus="{status}"><ClassifiedRecord>'
            f'<SimpleAllele AlleleID="{number}"><SequenceLocation Assembly="GRCh38" '
            f'Chr="1" positionVCF="{pos}" referenceAlleleVCF="{ref}" '
            f'alternateAlleleVCF="{alt}"/></SimpleAllele>'
            f'<Classifications><{classification} DateLastEvaluated="2020-01-01">'
            f"<Description>{term}</Description><ReviewStatus>criteria provided, single submitter"
            f"</ReviewStatus></{classification}></Classifications>"
            f'<ClinicalAssertion Accession="SCV90000000{number}" Version="1" '
            f'RecordStatus="{status}"><Submitter Name="Fabricated Lab"/>'
            f'<{classification} DateLastEvaluated="2020-01-01"><Description>{term}'
            f"</Description><ReviewStatus>criteria provided, single submitter</ReviewStatus>"
            f"</{classification}></ClinicalAssertion></ClassifiedRecord></VariationArchive>"
        )
    xml.write_text(
        '<ReleaseSet Dated="2026-07-01" ReleaseID="synthetic-m4">' + "".join(body) + "</ReleaseSet>"
    )
    ingest_clinvar_vcv(xml, tmp_path / "evidence")
    link_external_evidence(tmp_path / "m2", tmp_path / "evidence", tmp_path / "annotation")
    return tmp_path / "m2", tmp_path / "evidence", tmp_path / "annotation"


def _prioritize(tmp_path: Path) -> PrioritizationResult:
    m2, evidence, annotation = _runs(tmp_path)
    result = prioritize_clinical_variants(
        m2,
        evidence,
        annotation,
        tmp_path / "m4",
        ClinicalPrioritizationConfig(
            policy_path=POLICY, analysis_context=AnalysisContext.GERMLINE_CONSTITUTIONAL
        ),
    )
    return result


def test_priority_genotype_and_evidence_semantics(tmp_path: Path) -> None:
    result = _prioritize(tmp_path)
    by_pos = {profile.position: profile for profile in result.profiles}
    candidates = {
        by_pos[next(p for p, v in by_pos.items() if v.profile_id == c.profile_id)].position: c
        for c in result.candidates
    }
    assert candidates[10].priority_band == ReviewPriorityBand.REVIEW_FIRST
    assert candidates[20].priority_band == ReviewPriorityBand.NOT_ELIGIBLE
    assert by_pos[20].genotype_state == GenotypeEvidenceState.OBSERVED_REFERENCE_ONLY
    assert by_pos[30].genotype_state == GenotypeEvidenceState.NO_CANONICAL_CALLED_GENOTYPE
    assert candidates[30].priority_band == ReviewPriorityBand.NOT_ELIGIBLE
    assert by_pos[40].genotype_state == GenotypeEvidenceState.DISCORDANT_CALLED_ROWS_WITH_ALT
    assert candidates[40].eligibility == CandidateEligibility.DATA_CONFLICT
    assert candidates[40].priority_band == ReviewPriorityBand.REVIEW_NEXT
    assert len(by_pos[40].genotype_rows) == 2
    assert candidates[50].priority_band == ReviewPriorityBand.NOT_PRIORITIZED
    assert candidates[60].priority_band == ReviewPriorityBand.CONTEXT_ONLY
    assert all(age == 2373 for age in by_pos[10].evidence_age_days)
    assert set(by_pos[10].scv_assertion_ids).isdisjoint(by_pos[10].vcv_assertion_ids)
    assert any(r.reason_code == "data_conflict_alt_carriage_unresolved" for r in result.rationales)


def test_artifacts_policy_determinism_report_and_cli(tmp_path: Path) -> None:
    result = _prioritize(tmp_path)
    manifest = json.loads((tmp_path / "m4/manifest.json").read_text())
    assert result.policy_identity.file_sha256 == sha256(POLICY.read_bytes()).hexdigest()
    assert all(
        sha256((tmp_path / "m4" / name).read_bytes()).hexdigest() == digest
        for name, digest in manifest["artifacts"].items()
    )
    assert set(
        pl.read_parquet(tmp_path / "m4/candidate_assertion_links.parquet")["assertion_level"]
    )
    report = (tmp_path / "m4/prioritization_report.md").read_text()
    assert "not a negative genetic test" in report
    assert "No medical decisions" in report
    assert "SCV900000001.1" in report and "VCV900000001.1" in report
    m2, evidence, annotation = tmp_path / "m2", tmp_path / "evidence", tmp_path / "annotation"
    cli = CliRunner().invoke(
        app,
        [
            "prioritize",
            "clinical",
            "--normalization-run",
            str(m2),
            "--evidence-run",
            str(evidence),
            "--annotation-run",
            str(annotation),
            "--policy",
            str(POLICY),
            "--analysis-context",
            "germline_constitutional",
            "--output",
            str(tmp_path / "cli"),
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert "Pathogenic" not in cli.output and "SCV" not in cli.output
    second = CliRunner().invoke(
        app,
        [
            "prioritize",
            "clinical",
            "--normalization-run",
            str(m2),
            "--evidence-run",
            str(evidence),
            "--annotation-run",
            str(annotation),
            "--policy",
            str(POLICY),
            "--analysis-context",
            "germline_constitutional",
            "--output",
            str(tmp_path / "cli"),
        ],
    )
    assert second.exit_code == 2


def test_tampering_and_policy_validation_fail_atomically(tmp_path: Path) -> None:
    m2, evidence, annotation = _runs(tmp_path)
    link_path = annotation / "variant_evidence_links.parquet"
    original_links = link_path.read_bytes()
    link_path.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="checksum"):
        prioritize_clinical_variants(
            m2,
            evidence,
            annotation,
            tmp_path / "bad",
            ClinicalPrioritizationConfig(
                policy_path=POLICY, analysis_context=AnalysisContext.GERMLINE_CONSTITUTIONAL
            ),
        )
    assert not (tmp_path / "bad").exists()
    link_path.write_bytes(original_links)
    policy = json.loads(POLICY.read_text())
    policy["term_routes"]["benign_like"].append("Pathogenic")
    bad_policy = tmp_path / "bad-policy.json"
    bad_policy.write_text(json.dumps(policy))
    with pytest.raises(ValueError, match="invalid prioritization policy"):
        prioritize_clinical_variants(
            m2,
            evidence,
            annotation,
            tmp_path / "bad2",
            ClinicalPrioritizationConfig(
                policy_path=bad_policy, analysis_context=AnalysisContext.GERMLINE_CONSTITUTIONAL
            ),
        )
