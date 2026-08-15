import gzip
import io
import json
import subprocess
import threading
from hashlib import md5, sha256
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


def test_query_batch_does_not_retry_deterministically_rejected_output(tmp_path: Path) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())
    calls = 0

    def runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        output = Path(args[-1])
        output.write_text("malformed\n")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    with pytest.raises(resource_module._OutputValidationInterruption):  # noqa: SLF001
        resource_module._query_batch(  # noqa: SLF001
            tool,
            "a" * 64,
            (resource_module.DBSNP_URLS["GRCh37"],),
            "GRCh37",
            ("rs1",),
            tmp_path / "checkpoint",
            tmp_path / "work",
            runner,
            reporter,
            attempts=4,
            timeout_seconds=60,
            sleep=lambda _seconds: None,
        )

    assert calls == 1
    assert '"validation_category": "row-schema-invalid"' in reporter.log_path.read_text()
    assert "rs1" not in reporter.log_path.read_text()


def _common_identity(content: bytes = b"common-v1") -> resource_module.QueryResourceIdentity:
    return resource_module.QueryResourceIdentity(
        resource_module.QUERY_RESOURCE_IDENTITY_SCHEMA,
        "local-common-bigbed",
        "GRCh37",
        resource_module.DBSNP_BUILD,
        resource_module.DBSNP_COMMON_URLS["GRCh37"],
        sha256(content).hexdigest(),
        len(content),
    )


def test_common_checkpoint_identity_ignores_ephemeral_execution_path(tmp_path: Path) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())
    calls = 0

    def runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        Path(args[-1]).write_text("")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    for runtime in ("runtime-one", "runtime-two"):
        resource_module._query_batch(  # noqa: SLF001
            tool,
            "a" * 64,
            (f"/content/genome-evidence-resources-{runtime}/common/grch37/dbSnp155Common.bb",),
            "GRCh37",
            ("rs1",),
            tmp_path / "checkpoint",
            tmp_path / runtime,
            runner,
            reporter,
            attempts=1,
            timeout_seconds=60,
            sleep=lambda _seconds: None,
            allow_validation_indeterminate=True,
            indeterminate=[],
            resource_identity=_common_identity(),
        )

    assert calls == 1
    manifest = json.loads((tmp_path / "checkpoint/COMPLETED.json").read_text())
    assert manifest["resource_identity"] == _common_identity().as_dict()
    assert "/content" not in json.dumps(manifest)


def test_changed_common_content_invalidates_query_checkpoint(tmp_path: Path) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    calls = 0

    def runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        Path(args[-1]).write_text("")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    for content in (b"common-v1", b"common-v2"):
        resource_module._query_batch(  # noqa: SLF001
            tool,
            "a" * 64,
            ("/tmp/execution.bb",),
            "GRCh37",
            ("rs1",),
            tmp_path / "checkpoint",
            tmp_path / "work",
            runner,
            ProvisioningReporter(tmp_path / f"{calls}.jsonl", stream=io.StringIO()),
            attempts=1,
            timeout_seconds=60,
            sleep=lambda _seconds: None,
            resource_identity=_common_identity(content),
        )
    assert calls == 2


def test_legacy_common_checkpoint_is_narrowly_migrated(tmp_path: Path) -> None:
    directory = (
        tmp_path / "cache/downloads/normalization/v1/source-bound/dbsnp/common-query/"
        "grch37-common/batches/00000-synthetic"
    )
    directory.mkdir(parents=True)
    output = directory / "records.bed"
    output.write_text("")
    payload_digest = sha256(b"rs1\n").hexdigest()
    (directory / "COMPLETED.json").write_text(
        json.dumps(
            {
                "schema": resource_module.LEGACY_QUERY_CHECKPOINT_SCHEMA,
                "assembly": "GRCh37",
                "dbsnp_build": resource_module.DBSNP_BUILD,
                "canonical_url": (
                    "/content/genome-evidence-resources-old/common/grch37/dbSnp155Common.bb"
                ),
                "tool_sha256": "a" * 64,
                "identifiers_sha256": payload_digest,
                "identifier_count": 1,
                "record_count": 0,
                "sha256": sha256(b"").hexdigest(),
                "byte_size": 0,
            }
        )
    )
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())
    result = resource_module._valid_query_checkpoint(  # noqa: SLF001
        directory,
        identifiers=("rs1",),
        identifiers_sha256=payload_digest,
        assembly="GRCh37",
        resource_identity=_common_identity(),
        tool_sha256="a" * 64,
        allow_legacy_common_migration=True,
        reporter=reporter,
    )
    assert result == output
    assert json.loads((directory / "COMPLETED.json").read_text())["schema"].endswith("/v2")
    assert "dbsnp.batch.migrated" in reporter.log_path.read_text()


def test_common_indeterminate_leaf_persists_but_strict_full_does_not(tmp_path: Path) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    calls = 0

    def malformed(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        Path(args[-1]).write_text("malformed\n")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    for _ in range(2):
        leaves: list[tuple[str, ...]] = []
        assert (
            resource_module._query_batch(  # noqa: SLF001
                tool,
                "a" * 64,
                ("/tmp/common.bb",),
                "GRCh37",
                ("rs1",),
                tmp_path / "checkpoint",
                tmp_path / "work",
                malformed,
                ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO()),
                attempts=1,
                timeout_seconds=60,
                sleep=lambda _seconds: None,
                allow_validation_indeterminate=True,
                indeterminate=leaves,
                resource_identity=_common_identity(),
            )
            is None
        )
        assert leaves == [("rs1",)]
    assert calls == 1
    with pytest.raises(resource_module._OutputValidationInterruption):  # noqa: SLF001
        resource_module._query_batch(  # noqa: SLF001
            tool,
            "a" * 64,
            (resource_module.DBSNP_URLS["GRCh37"],),
            "GRCh37",
            ("rs1",),
            tmp_path / "strict",
            tmp_path / "work",
            malformed,
            ProvisioningReporter(tmp_path / "strict.jsonl", stream=io.StringIO()),
            attempts=1,
            timeout_seconds=60,
            sleep=lambda _seconds: None,
        )
    assert not (tmp_path / "strict/INDETERMINATE.json").exists()


def test_common_validation_failure_is_localized_to_leaf_and_missing_ids(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(resource_module, "_MIN_SPLIT_BATCH_SIZE", 2)
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())
    identifiers = tuple(f"rs{index}" for index in range(1, 7))
    full_requests: list[tuple[str, ...]] = []

    def download(_url: str, destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(b"synthetic common BigBed")

    def runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        requested = tuple(Path(args[-2]).read_text().splitlines())
        output = Path(args[-1])
        if args[2] == resource_module.DBSNP_URLS["GRCh37"]:
            full_requests.append(requested)
            output.write_text("")
        elif "dbSnp155Common.bb" in args[2]:
            if "rs3" in requested:
                output.write_text("malformed\n")
            else:
                returned = [value for value in requested if value in {"rs1", "rs5"}]
                output.write_text(
                    "".join(f"chr1\t1\t2\t{value}\tC\t1\tG,\t0\n" for value in returned)
                )
        else:
            raise AssertionError("unexpected query endpoint")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    _, provenance = resource_module._query_common_first(  # noqa: SLF001
        tool,
        "a" * 64,
        "GRCh37",
        identifiers,
        tmp_path / "checkpoints",
        tmp_path / "work",
        runner,
        reporter,
        download=download,
        batch_size=4,
        attempts=4,
        timeout_seconds=60,
        workers=1,
        download_segments=1,
        sleep=lambda _seconds: None,
        label="grch37",
    )

    assert full_requests == [("rs2", "rs3", "rs4", "rs6")]
    assert provenance["common_hit_count"] == 2
    assert provenance["common_missing_count"] == 2
    assert provenance["common_indeterminate_count"] == 2
    assert provenance["full_fallback_requested_count"] == 4
    assert provenance["completed_full_fallback_count"] == 4
    assert "rs3" not in reporter.log_path.read_text()


def test_common_fallback_plan_stays_proportional_at_real_world_scale() -> None:
    identifiers = tuple(f"rs{index}" for index in range(950_000))
    indeterminate = identifiers[320_000:320_156]
    missing = identifiers[-1_000:]
    indeterminate_set = set(indeterminate)
    missing_set = set(missing)
    unresolved = indeterminate_set | missing_set
    covered = tuple(value for value in identifiers if value not in indeterminate_set)
    returned = tuple(value for value in covered if value not in missing_set)
    result = resource_module._CommonQueryResult(  # noqa: SLF001
        output=None,
        validated_outputs=(),
        covered_identifiers=covered,
        returned_identifiers=returned,
        missing_identifiers=missing,
        indeterminate_identifiers=indeterminate,
    )

    fallback = resource_module._common_fallback_identifiers(identifiers, result)  # noqa: SLF001

    assert len(fallback) == 1_156
    assert set(fallback) == unresolved
    assert len(fallback) < len(identifiers) // 100


def test_transient_nonzero_query_failure_retries_then_completes(tmp_path: Path) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    calls = 0

    def runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        if calls == 1:
            return subprocess.CompletedProcess(args, 255, stdout="", stderr="temporary failure")
        Path(args[-1]).write_text("")
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
        ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO()),
        attempts=2,
        timeout_seconds=60,
        sleep=lambda _seconds: None,
    )

    assert output is not None
    assert calls == 2


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


def test_full_fallback_batches_overlap_with_isolated_worker_caches(tmp_path: Path) -> None:
    tool = tmp_path / "bigBedNamedItems"
    tool.write_bytes(b"synthetic")
    identifiers = tuple(f"rs{index}" for index in range(500))
    barrier = threading.Barrier(2)
    lock = threading.Lock()
    active = 0
    maximum = 0
    tmpdirs: set[str] = set()

    def runner(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        nonlocal active, maximum
        environment = kwargs["env"]
        assert isinstance(environment, dict)
        tmpdirs.add(str(environment["TMPDIR"]))
        with lock:
            active += 1
            maximum = max(maximum, active)
        barrier.wait(timeout=2)
        Path(args[-1]).write_text("")
        with lock:
            active -= 1
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    output = resource_module._query_bigbed_in_batches(  # noqa: SLF001
        tool,
        "a" * 64,
        (resource_module.DBSNP_URLS["GRCh37"],),
        "GRCh37",
        identifiers,
        tmp_path / "checkpoint",
        tmp_path / "work",
        runner,
        ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO()),
        batch_size=250,
        attempts=1,
        timeout_seconds=60,
        sleep=lambda _seconds: None,
        label="synthetic",
        workers=2,
    )

    assert output.read_text() == ""
    assert maximum == 2
    assert len(tmpdirs) == 2
    assert all((Path(directory) / "udcCache").is_dir() for directory in tmpdirs)
