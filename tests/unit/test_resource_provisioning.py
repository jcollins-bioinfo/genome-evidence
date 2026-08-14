import gzip
import io
import json
import subprocess
from hashlib import md5
from pathlib import Path

import pytest

import genome_evidence.workspace.resources as resource_module
from genome_evidence.workspace.provisioning_progress import ProvisioningReporter
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


def test_fasta_provisioning_repairs_a_torn_completion_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = gzip.compress(b">chr1\nACGT\n")
    monkeypatch.setattr(
        resource_module,
        "FASTA_UPSTREAM_MD5",
        md5(payload, usedforsecurity=False).hexdigest(),
    )
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())

    def download(_url: str, destination: Path) -> None:
        destination.write_bytes(payload)

    fasta, index, manifest = resource_module._prepare_fasta(  # noqa: SLF001
        tmp_path / "workspace",
        tmp_path / "checkpoint",
        download,
        reporter,
        sleep=lambda _seconds: None,
    )
    completed = fasta.parent / "COMPLETED.json"
    completed.write_text("{", encoding="utf-8")

    repaired_fasta, repaired_index, repaired_manifest = resource_module._prepare_fasta(  # noqa: SLF001
        tmp_path / "workspace",
        tmp_path / "checkpoint",
        download,
        reporter,
        sleep=lambda _seconds: None,
    )

    assert repaired_fasta == fasta
    assert repaired_index == index
    assert repaired_manifest == manifest
    assert json.loads(completed.read_text()) == manifest


def test_kent_tool_probe_rejects_an_unrelated_executable(tmp_path: Path) -> None:
    root = tmp_path / "workspace"
    cached = root / "cache/tools/ucsc/kent-v479/bigBedNamedItems"
    cached.parent.mkdir(parents=True)
    cached.write_bytes(b"not the expected tool")
    work = tmp_path / "work"
    work.mkdir()

    def unused_download(_url: str, _destination: Path) -> None:
        raise AssertionError("existing cache must not be downloaded again")

    def wrong_runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(args, 1, stdout="", stderr="unrelated usage\n")

    with pytest.raises(ValueError, match="expected bigBedNamedItems CLI"):
        resource_module._prepare_kent_tool(  # noqa: SLF001
            root,
            work,
            unused_download,
            wrong_runner,
            ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO()),
            sleep=lambda _seconds: None,
        )


def test_query_batch_retries_and_rejects_malformed_zero_exit_output(tmp_path: Path) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())
    calls = 0

    def runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        output = Path(args[-1])
        output.write_text("malformed\n" if calls == 1 else "chr1\t1\t2\trs1\tC\t1\tG,\t0\n")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    output = resource_module._query_batch(  # noqa: SLF001
        tool,
        "a" * 64,
        (resource_module.DBSNP_URLS["GRCh37"],),
        "GRCh37",
        ("rs1",),
        tmp_path / "checkpoint",
        tmp_path / "work",
        runner,
        reporter,
        attempts=2,
        timeout_seconds=60,
        sleep=lambda _seconds: None,
    )

    assert calls == 2
    assert "rs1" in output.read_text()
    assert "rs1" not in reporter.log_path.read_text()


def test_query_batch_resumes_a_persisted_split_without_repeating_parent(
    tmp_path: Path,
) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())
    identifiers = tuple(f"rs{index}" for index in range(1, 301))
    queried_counts: list[int] = []

    def failing_runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        queried_counts.append(len(Path(args[-2]).read_text().splitlines()))
        return subprocess.CompletedProcess(
            args,
            255,
            stdout="",
            stderr="remote failure for rsSensitive",
        )

    arguments = (
        tool,
        "a" * 64,
        (resource_module.DBSNP_URLS["GRCh37"],),
        "GRCh37",
        identifiers,
        tmp_path / "checkpoint",
        tmp_path / "work",
        failing_runner,
        reporter,
    )
    with pytest.raises(resource_module._OperationalInterruption):  # noqa: SLF001
        resource_module._query_batch(  # noqa: SLF001
            *arguments,
            attempts=1,
            timeout_seconds=60,
            sleep=lambda _seconds: None,
        )
    assert queried_counts == [300, 150]

    queried_counts.clear()
    with pytest.raises(resource_module._OperationalInterruption):  # noqa: SLF001
        resource_module._query_batch(  # noqa: SLF001
            *arguments,
            attempts=1,
            timeout_seconds=60,
            sleep=lambda _seconds: None,
        )

    assert queried_counts == [150]
    assert "split.resume" in reporter.log_path.read_text()
    assert "rsSensitive" not in reporter.log_path.read_text()


def test_bigbed_value_model_is_immutable() -> None:
    value = BigBedVariant("1", 1, "rs1", "A", ("G",))
    with pytest.raises(AttributeError):
        value.position = 2  # type: ignore[misc]
