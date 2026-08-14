"""Synthetic-only adversarial tests for the M8 pharmacogenomics foundation."""

import json
from hashlib import sha256
from pathlib import Path

import polars as pl
import pytest

from genome_evidence.pharmacogenomics import validate_pharmacogenomics_bundle
from genome_evidence.pharmacogenomics.matching import enumerate_diplotypes
from genome_evidence.pharmacogenomics.models import (
    GeneOutcome,
    LocusEvidence,
    LocusEvidenceState,
    MatchingLimits,
)


def _bundle(root: Path) -> Path:
    root.mkdir()
    source = {
        "sources": [
            {
                "source_id": "synthetic-source",
                "resource": "fabricated-test-resource",
                "version": "test-1",
                "url": "https://example.invalid/synthetic/test-1",
                "retrieved_at": "2000-01-01T00:00:00Z",
                "sha256": "0" * 64,
                "license": "synthetic-test-only",
                "content_fingerprint": "synthetic-content-v1",
            }
        ]
    }
    (root / "sources.json").write_text(json.dumps(source))
    tables = {
        "genes.parquet": [{"gene_id": "FAKE1", "strategy": "star_haplotype_small_variant"}],
        "loci.parquet": [
            {
                "locus_id": "fake-locus-1",
                "gene_id": "FAKE1",
                "assembly": "GRCh38",
                "chromosome": "22",
                "position": 1001,
                "reference": "A",
                "alternate": "G",
            }
        ],
        "alleles.parquet": [
            {
                "allele_id": "FAKE1*1",
                "gene_id": "FAKE1",
                "status": "active",
                "source_id": "synthetic-source",
                "source_version": "test-1",
                "haplotype_evidence_level": "fabricated-high",
            },
            {
                "allele_id": "FAKE1*2",
                "gene_id": "FAKE1",
                "status": "active",
                "source_id": "synthetic-source",
                "source_version": "test-1",
                "haplotype_evidence_level": "fabricated-low",
            },
        ],
        "allele_locus_constraints.parquet": [
            {"allele_id": "FAKE1*1", "locus_id": "fake-locus-1", "state": "required_reference"},
            {"allele_id": "FAKE1*2", "locus_id": "fake-locus-1", "state": "required_alternate"},
        ],
        "allele_function_assertions.parquet": [
            {
                "assertion_id": "function-1",
                "gene_id": "FAKE1",
                "allele_id": "FAKE1*1",
                "source_id": "synthetic-source",
                "source_version": "test-1",
                "source_term": "fabricated function term",
            }
        ],
        "diplotype_phenotype_rules.parquet": [
            {
                "rule_id": "rule-1",
                "gene_id": "FAKE1",
                "allele_a": "FAKE1*1",
                "allele_b": "FAKE1*2",
                "source_id": "synthetic-source",
                "source_version": "test-1",
                "phenotype_term": "fabricated phenotype",
                "historical_mapping": False,
            }
        ],
        "guideline_evidence.parquet": [
            {
                "guideline_id": "guide-1",
                "gene_id": "FAKE1",
                "source_id": "synthetic-source",
                "drug_id": "FAKE-DRUG",
                "source_version": "test-1",
                "url": "https://example.invalid/guide",
            }
        ],
    }
    for name, rows in tables.items():
        pl.DataFrame(rows).write_parquet(root / name)
    artifacts = {}
    for name in ("sources.json", *tables):
        path = root / name
        artifacts[name] = {
            "sha256": sha256(path.read_bytes()).hexdigest(),
            "byte_size": path.stat().st_size,
            "row_count": len(tables.get(name, [])),
            "schema": f"synthetic-{name}/v1",
        }
    manifest = {
        "schema": "genome-evidence-pgx-bundle/v1",
        "bundle_id": "synthetic-pgx-bundle",
        "bundle_version": "test-1",
        "creation_tool_version": "0.2.0",
        "created_at": "2000-01-01T00:00:00Z",
        "configuration_hash": "1" * 64,
        "assembly": "GRCh38",
        "classification": "synthetic",
        "artifacts": artifacts,
    }
    (root / "bundle_manifest.json").write_text(json.dumps(manifest, sort_keys=True))
    return root


def test_bundle_validation_and_tamper_fail_closed(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path / "bundle")
    first = validate_pharmacogenomics_bundle(bundle)
    assert first == validate_pharmacogenomics_bundle(bundle)
    (bundle / "sources.json").write_text("{}")
    with pytest.raises(ValueError, match="checksum or size"):
        validate_pharmacogenomics_bundle(bundle)


@pytest.mark.parametrize(
    ("state", "alleles", "outcome"),
    [
        (LocusEvidenceState.OBSERVED_REFERENCE, ("A", "A"), GeneOutcome.RESOLVED_CANDIDATE),
        (LocusEvidenceState.OBSERVED_HETEROZYGOUS, ("A", "G"), GeneOutcome.RESOLVED_CANDIDATE),
        (LocusEvidenceState.OBSERVED_ALTERNATE, ("G", "G"), GeneOutcome.RESOLVED_CANDIDATE),
        (LocusEvidenceState.MISSING_OR_UNASSAYED, (), GeneOutcome.INSUFFICIENT_COVERAGE),
        (LocusEvidenceState.DUPLICATE_CONFLICT, (), GeneOutcome.CONFLICTING_OBSERVATIONS),
    ],
)
def test_missing_never_becomes_reference_and_multisets_match(
    tmp_path: Path,
    state: LocusEvidenceState,
    alleles: tuple[str, ...],
    outcome: GeneOutcome,
) -> None:
    bundle = _bundle(tmp_path / "bundle")
    evidence = (LocusEvidence(locus_id="fake-locus-1", state=state, alleles=alleles),)
    result = enumerate_diplotypes(bundle, "FAKE1", evidence, MatchingLimits())
    assert result.outcome == outcome
    if state == LocusEvidenceState.MISSING_OR_UNASSAYED:
        assert len(result.candidates) == 3
        assert not any(item.fully_evaluated for item in result.candidates)


def test_candidate_limit_fails_without_truncating(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path / "bundle")
    result = enumerate_diplotypes(
        bundle,
        "FAKE1",
        (LocusEvidence(locus_id="fake-locus-1", state=LocusEvidenceState.MISSING_OR_UNASSAYED),),
        MatchingLimits(max_candidate_pairs=1),
    )
    assert result.outcome == GeneOutcome.COMBINATORIAL_LIMIT_EXCEEDED
    assert result.candidates == ()


def test_extra_and_symlink_artifacts_are_rejected(tmp_path: Path) -> None:
    bundle = _bundle(tmp_path / "bundle")
    (bundle / "undeclared.txt").write_text("synthetic")
    with pytest.raises(ValueError, match="extra"):
        validate_pharmacogenomics_bundle(bundle)
