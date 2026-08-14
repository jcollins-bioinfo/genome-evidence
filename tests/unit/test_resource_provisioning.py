import gzip
from pathlib import Path

import pytest

from genome_evidence.workspace.resources import (
    BigBedVariant,
    SourceMarker,
    build_fasta_index,
    build_marker_resources,
    parse_bigbed_variants,
    read_source_markers,
)


def test_source_reader_retains_only_marker_identity_and_coordinates(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("# genome build: GRCh37\nrs1\t1\t2\tCG\ni1\tX\t9\t--\n")

    assembly, markers = read_source_markers(source)

    assert assembly == "GRCh37"
    assert markers == (SourceMarker("rs1", "1", 2), SourceMarker("i1", "X", 9))
    assert all(not hasattr(marker, "genotype") for marker in markers)


@pytest.mark.parametrize(("vendor_build", "expected"), [("37", "GRCh37"), ("38", "GRCh38")])
def test_source_reader_recognizes_standard_23andme_build_header(
    tmp_path: Path, vendor_build: str, expected: str
) -> None:
    source = tmp_path / "source.txt"
    source.write_text(
        f"# We are using reference human assembly build {vendor_build}.\nrs1\t1\t2\tCG\n"
    )

    assert read_source_markers(source)[0] == expected


def test_source_reader_requires_a_verified_build(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("rs1\t1\t2\tCG\n")

    with pytest.raises(ValueError, match="source build is missing"):
        read_source_markers(source)
    assert read_source_markers(source, "hg38")[0] == "GRCh38"


def test_bigbed_parser_and_resource_builder_are_exact_and_fail_closed() -> None:
    source_rows = parse_bigbed_variants(
        "chr1\t1\t2\trs1\tC\t1\tG,\t0\n"
        "chr1\t6\t7\trs2\tA\t2\tC,T,\t0\n"
        "chr1\t8\t9\trsWrongPosition\tG\t1\tA,\t0\n"
    )
    target_rows = parse_bigbed_variants(
        "chr1\t3\t4\trs1\tC\t1\tG,\t0\nchr2\t10\t11\trs2\tA\t2\tC,T,\t0\n"
    )

    definitions, mappings, stats = build_marker_resources(
        (
            SourceMarker("rs1", "1", 2),
            SourceMarker("rs2", "1", 7),
            SourceMarker("rsWrongPosition", "1", 10),
            SourceMarker("internal", "1", 20),
        ),
        "GRCh37",
        source_rows,
        target_rows,
    )

    assert [row["marker_id"] for row in definitions] == ["rs1", "rs2", "rs2"]
    assert all(row["orientation"] == "none" for row in definitions)
    assert all(row["orientation_authoritative"] is True for row in definitions)
    assert mappings == {"1:1": [["1", 3]], "1:6": [["2", 10]]}
    assert stats == {
        "source_marker_count": 4,
        "rsid_marker_count": 3,
        "defined_marker_count": 2,
        "definition_count": 3,
        "cross_build_mapped_marker_count": 2,
        "exact_source_placement_count": 2,
    }


def test_fasta_index_supports_bounded_random_access(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_bytes(b">chr1 description\nACGT\nAC\n>chrM\nTTTT\n")
    index = tmp_path / "reference.fa.fai"

    build_fasta_index(fasta, index)

    assert index.read_text() == "chr1\t6\t18\t4\t5\nchrM\t4\t32\t4\t5\n"


def test_fasta_index_rejects_an_interior_short_line(tmp_path: Path) -> None:
    fasta = tmp_path / "bad.fa"
    with gzip.open(tmp_path / "unused.gz", "wb"):
        pass
    fasta.write_bytes(b">chr1\nAAAA\nAA\nAAAA\n")

    with pytest.raises(ValueError, match="irregular non-terminal"):
        build_fasta_index(fasta, tmp_path / "bad.fa.fai")


def test_bigbed_value_model_is_immutable() -> None:
    value = BigBedVariant("1", 1, "rs1", "A", ("G",))
    with pytest.raises(AttributeError):
        value.position = 2  # type: ignore[misc]
