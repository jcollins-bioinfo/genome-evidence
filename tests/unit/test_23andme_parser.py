from hashlib import sha256
from pathlib import Path

import pytest
from pydantic import ValidationError

from genome_evidence.ingest import Ingest23andMeConfig, ParseMode, ingest_23andme
from genome_evidence.ingest.errors import GenotypeParseError
from genome_evidence.qc.models import BuildProvenance, CallState, LexicalGenotypeCategory

FIXTURES = Path(__file__).parents[1] / "fixtures" / "23andme"


def ingest(name: str, tmp_path: Path, **config: object):
    return ingest_23andme(FIXTURES / name, tmp_path / "run", Ingest23andMeConfig(**config))


def test_normal_source_is_preserved_and_accounted(tmp_path: Path) -> None:
    result = ingest("normal.txt", tmp_path)
    first = result.observations[0]
    assert (
        first.source_marker_id,
        first.source_chromosome,
        first.source_position,
        first.raw_genotype,
    ) == ("rsSynthetic1", "1", 101, "AG")
    assert result.source_metadata.comment_lines[0].startswith("# SYNTHETIC")
    assert result.source_metadata.declared_build == "GRCh37"
    assert (
        result.qc_summary.physical_line_count
        == result.qc_summary.comment_line_count
        + result.qc_summary.blank_line_count
        + result.qc_summary.data_line_count
    )
    assert (
        result.qc_summary.eligible_parsed_record_count
        == result.qc_summary.called_record_count + result.qc_summary.no_call_record_count
    )
    assert result.qc_summary.call_rate == pytest.approx(2 / 3)
    chromosome_one = result.qc_summary.chromosome_summaries[0]
    assert (
        chromosome_one.record_count,
        chromosome_one.called_count,
        chromosome_one.min_source_position,
        chromosome_one.max_source_position,
    ) == (2, 2, 101, 202)
    assert not hasattr(first, "reference") and not hasattr(first, "alternate")


def test_sha256_is_exact_and_deterministic(tmp_path: Path) -> None:
    path = FIXTURES / "normal.txt"
    expected = sha256(path.read_bytes()).hexdigest()
    assert ingest("normal.txt", tmp_path).qc_summary.source_sha256 == expected


@pytest.mark.parametrize(
    "name,marker", [("crlf.txt", "syntheticCRLF"), ("bom.txt", "syntheticBOM")]
)
def test_transport_encodings(name: str, marker: str, tmp_path: Path) -> None:
    assert ingest(name, tmp_path).observations[0].source_marker_id == marker


def test_no_call_is_not_a_called_reference(tmp_path: Path) -> None:
    records = ingest("normal.txt", tmp_path).observations
    no_call = records[-1]
    assert no_call.call_status == CallState.NO_CALL
    assert no_call.raw_genotype == "--"
    assert records[1].call_status == CallState.CALLED
    assert len(records) == 3  # no unobserved locus is synthesized


@pytest.mark.parametrize(
    "name,metric",
    [
        ("duplicate_marker.txt", "duplicate_marker_id_count"),
        ("duplicate_coordinate.txt", "duplicate_coordinate_count"),
    ],
)
def test_duplicates_are_preserved_and_reported(name: str, metric: str, tmp_path: Path) -> None:
    result = ingest(name, tmp_path)
    assert len(result.observations) == 2
    assert getattr(result.qc_summary, metric) == 1


def test_conflicting_duplicate_is_reported(tmp_path: Path) -> None:
    result = ingest("conflicting_duplicate.txt", tmp_path)
    assert result.qc_summary.conflicting_duplicate_marker_count == 1
    assert any(f.code == "CONFLICTING_DUPLICATE_MARKER" for f in result.findings)


def test_strict_rejects_malformed_without_genotype_dump(tmp_path: Path) -> None:
    with pytest.raises(GenotypeParseError, match="line 2: Source position") as captured:
        ingest("malformed_position.txt", tmp_path)
    assert "syntheticBad" not in str(captured.value)


def test_lenient_continues_and_accounts_for_malformed(tmp_path: Path) -> None:
    result = ingest("malformed_position.txt", tmp_path, mode=ParseMode.LENIENT)
    assert [r.source_marker_id for r in result.observations] == ["syntheticGood"]
    assert result.qc_summary.malformed_record_count == 1
    assert result.qc_summary.invalid_position_count == 1
    assert result.qc_summary.data_line_count == 2
    assert result.malformed_records[0].source_line_number == 2


@pytest.mark.parametrize("name", ["malformed_columns.txt", "malformed_position.txt"])
def test_structural_errors_fail_in_strict_mode(name: str, tmp_path: Path) -> None:
    with pytest.raises(GenotypeParseError):
        ingest(name, tmp_path)


def test_unknown_values_are_preserved_and_lexically_described(tmp_path: Path) -> None:
    unknown = ingest("unknown_chromosome.txt", tmp_path).observations[0]
    assert unknown.source_chromosome == "vendor_contig"
    assert unknown.raw_genotype == "AC"
    result = ingest("unusual_genotype.txt", tmp_path / "other")
    assert [r.raw_genotype for r in result.observations] == ["A", "DEL"]
    assert result.observations[0].lexical_category == LexicalGenotypeCategory.SINGLE_ALLELE_TOKEN
    assert result.observations[1].lexical_category == LexicalGenotypeCategory.OTHER_TOKEN


def test_non_rs_identifier_is_unchanged(tmp_path: Path) -> None:
    assert (
        ingest("non_rs_marker.txt", tmp_path).observations[0].source_marker_id
        == "vendor-internal-id"
    )


def test_build_is_not_guessed_and_override_is_tagged(tmp_path: Path) -> None:
    unknown = ingest("missing_build.txt", tmp_path)
    assert unknown.qc_summary.declared_or_resolved_build == "UNKNOWN"
    assert unknown.qc_summary.build_provenance == BuildProvenance.UNKNOWN
    overridden = ingest("missing_build.txt", tmp_path / "override", genome_build_override="GRCh38")
    assert overridden.source_metadata.declared_build is None
    assert overridden.source_metadata.resolved_build == "GRCh38"
    assert overridden.source_metadata.build_provenance == BuildProvenance.USER_OVERRIDE


def test_out_of_order_is_reported_without_reordering(tmp_path: Path) -> None:
    result = ingest("out_of_order.txt", tmp_path)
    assert [r.source_position for r in result.observations] == [100, 90]
    assert result.qc_summary.out_of_order_record_count == 1


def test_empty_source_has_defined_call_rate(tmp_path: Path) -> None:
    result = ingest("empty.txt", tmp_path)
    assert result.observations == ()
    assert result.qc_summary.call_rate is None


def test_source_observation_is_immutable(tmp_path: Path) -> None:
    record = ingest("normal.txt", tmp_path).observations[0]
    with pytest.raises(ValidationError):
        record.raw_genotype = "TT"
