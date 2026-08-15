import http.client
import io
import json
import stat
import urllib.error
import urllib.request
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pytest

import genome_evidence.workspace.provisioning_progress as progress_module
from genome_evidence.workspace.provisioning_progress import (
    ProvisioningReporter,
    WorkflowProgress,
    progress_bar,
    resumable_download,
)


class _HTTPResponse:
    def __init__(self, payload: bytes, *, status: int, headers: dict[str, str]) -> None:
        self._body = io.BytesIO(payload)
        self.status = status
        self.headers = headers

    def read(self, size: int = -1) -> bytes:
        return self._body.read(size)

    def getcode(self) -> int:
        return self.status

    def __enter__(self) -> "_HTTPResponse":
        return self

    def __exit__(self, *_args: object) -> None:
        return None


class _InterruptedHTTPResponse(_HTTPResponse):
    def __init__(self) -> None:
        super().__init__(
            b"abc",
            status=200,
            headers={"Content-Length": "6", "ETag": '"reference-v1"'},
        )
        self._completed_first_read = False

    def read(self, size: int = -1) -> bytes:
        if not self._completed_first_read:
            self._completed_first_read = True
            return super().read(size)
        raise http.client.IncompleteRead(b"", 3)


def _events(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def _clock(values: Iterator[float]) -> Callable[[], float]:
    return lambda: next(values)


def test_progress_bars_have_stable_width_and_workflow_is_monotonic() -> None:
    for style in ("unicode", "ascii"):
        assert len(progress_bar(0.0, width=17, style=style)) == 17
        assert len(progress_bar(0.537, width=17, style=style)) == 17
        assert len(progress_bar(1.0, width=17, style=style)) == 17
    workflow = WorkflowProgress()
    values = [
        workflow.update("preflight", 0.5),
        workflow.complete("grch37_fallback"),
        workflow.update("preflight", 0.1),
        workflow.complete("preflight"),
    ]
    assert values == sorted(values)


def test_dashboard_failure_is_isolated_and_jsonl_matches_numeric_state(tmp_path: Path) -> None:
    calls = 0

    def broken_display(_value: str) -> None:
        nonlocal calls
        calls += 1
        raise RuntimeError("synthetic display failure")

    reporter = ProvisioningReporter(
        tmp_path / "events.jsonl",
        stream=io.StringIO(),
        display_sink=broken_display,
        progress_interval_seconds=0,
    )
    reporter.workflow_progress("preflight", 0.5, force=True, resumed_batches=2)
    reporter.workflow_progress("preflight", 1.0, force=True, resumed_batches=2)
    records = [row for row in _events(reporter.log_path) if row["event"] == "workflow.progress"]
    assert calls == 1
    assert records[-1]["stage_fraction"] == 1.0
    assert records[-1]["overall_fraction"] == pytest.approx(0.03)
    assert records[-1]["resumed_batches"] == 2


def test_reporter_mirrors_aggregate_progress_to_stdout_and_private_jsonl(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "logs" / "00b.jsonl"
    stdout = io.StringIO()
    reporter = ProvisioningReporter(
        log_path,
        stream=stdout,
        clock=_clock(iter((12.0,))),
        progress_interval_seconds=0,
    )

    reporter.info(
        "source.validated",
        "Validated imported source metadata.",
        assembly="GRCh37",
        marker_count=2,
    )
    reporter.progress(
        "download.reference",
        "hg38.fa.gz",
        4,
        8,
        unit="bytes",
        started_at=10.0,
        force=True,
    )
    reporter.success(
        "resource.complete",
        "Published checksummed normalization resources.",
        defined_marker_count=1,
        unresolved_marker_count=1,
    )

    console = stdout.getvalue()
    records = _events(log_path)
    assert "INFO  source.validated" in console
    assert "PROG  download.reference" in console
    assert "DONE  resource.complete" in console
    assert "50.0%" in console
    assert [record["event"] for record in records] == [
        "source.validated",
        "download.reference",
        "resource.complete",
    ]
    assert records[1]["completed"] == 4
    assert records[1]["total"] == 8
    assert records[1]["rate_per_second"] == pytest.approx(2.0)
    assert records[1]["eta_seconds"] == pytest.approx(2.0)
    assert stat.S_IMODE(log_path.stat().st_mode) == 0o600

    combined = console + log_path.read_text(encoding="utf-8")
    assert "rs123" not in combined
    assert "genotype" not in combined.lower()
    assert "\tAA" not in combined


def test_reporter_keeps_stdout_progress_when_the_log_sink_is_unavailable(
    tmp_path: Path,
) -> None:
    blocked_parent = tmp_path / "not-a-directory"
    blocked_parent.write_text("blocked")
    stdout = io.StringIO()
    reporter = ProvisioningReporter(blocked_parent / "events.jsonl", stream=stdout)

    reporter.info("source.validated", "Validated aggregate source metadata.")

    console = stdout.getvalue()
    assert "log.unavailable" in console
    assert "source.validated" in console


def test_resumable_download_continues_after_early_eof_with_http_range(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "cache" / "reference.bin"
    log_path = tmp_path / "events.jsonl"
    reporter = ProvisioningReporter(
        log_path,
        stream=io.StringIO(),
        progress_interval_seconds=0,
    )
    responses = iter(
        (
            _HTTPResponse(
                b"abc",
                status=200,
                headers={"Content-Length": "6", "ETag": '"reference-v1"'},
            ),
            _HTTPResponse(
                b"def",
                status=206,
                headers={
                    "Content-Length": "3",
                    "Content-Range": "bytes 3-5/6",
                    "ETag": '"reference-v1"',
                },
            ),
        )
    )
    requested_ranges: list[str | None] = []
    requested_validators: list[str | None] = []

    def fake_urlopen(request: urllib.request.Request, *, timeout: int) -> _HTTPResponse:
        assert timeout == 120
        requested_ranges.append(request.get_header("Range"))
        requested_validators.append(request.get_header("If-range"))
        return next(responses)

    monkeypatch.setattr(progress_module.urllib.request, "urlopen", fake_urlopen)

    result = resumable_download(
        "https://example.test/reference.bin",
        destination,
        reporter,
        attempts=2,
        sleep=lambda _seconds: None,
        chunk_bytes=2,
    )

    assert result == destination
    assert destination.read_bytes() == b"abcdef"
    assert not destination.with_name("reference.bin.part").exists()
    assert requested_ranges == [None, "bytes=3-"]
    assert requested_validators == [None, '"reference-v1"']
    records = _events(log_path)
    assert any(record["event"] == "download.retry" for record in records)
    second_attempt = [record for record in records if record["event"] == "download.attempt"][1]
    assert second_attempt["resume_offset"] == 3


def test_resumable_download_retries_an_explicit_incomplete_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "cache" / "reference.bin"
    reporter = ProvisioningReporter(
        tmp_path / "events.jsonl",
        stream=io.StringIO(),
        progress_interval_seconds=0,
    )
    responses = iter(
        (
            _InterruptedHTTPResponse(),
            _HTTPResponse(
                b"def",
                status=206,
                headers={
                    "Content-Length": "3",
                    "Content-Range": "bytes 3-5/6",
                    "ETag": '"reference-v1"',
                },
            ),
        )
    )

    def fake_urlopen(request: urllib.request.Request, *, timeout: int) -> _HTTPResponse:
        assert timeout == 120
        return next(responses)

    monkeypatch.setattr(progress_module.urllib.request, "urlopen", fake_urlopen)

    resumable_download(
        "https://example.test/reference.bin",
        destination,
        reporter,
        attempts=2,
        sleep=lambda _seconds: None,
    )

    assert destination.read_bytes() == b"abcdef"


def test_resumable_download_restarts_safely_when_server_ignores_range(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "cache" / "reference.bin"
    destination.parent.mkdir(parents=True)
    destination.with_name("reference.bin.part").write_bytes(b"stale-prefix")
    destination.with_name("reference.bin.part.json").write_text(
        json.dumps(
            {
                "url": "https://example.test/reference.bin",
                "etag": '"stale-v1"',
                "last_modified": None,
                "total_bytes": 17,
            }
        )
    )
    log_path = tmp_path / "events.jsonl"
    stdout = io.StringIO()
    reporter = ProvisioningReporter(log_path, stream=stdout, progress_interval_seconds=0)
    requested_ranges: list[str | None] = []

    def fake_urlopen(request: urllib.request.Request, *, timeout: int) -> _HTTPResponse:
        assert timeout == 120
        requested_ranges.append(request.get_header("Range"))
        return _HTTPResponse(b"fresh", status=200, headers={"Content-Length": "5"})

    monkeypatch.setattr(progress_module.urllib.request, "urlopen", fake_urlopen)

    resumable_download(
        "https://example.test/reference.bin",
        destination,
        reporter,
        attempts=1,
        sleep=lambda _seconds: None,
        chunk_bytes=2,
    )

    assert requested_ranges == ["bytes=12-"]
    assert destination.read_bytes() == b"fresh"
    assert "download.restart" in stdout.getvalue()
    restart = next(record for record in _events(log_path) if record["event"] == "download.restart")
    assert restart["http_status"] == 200


def test_resumable_download_accepts_416_only_for_an_exact_complete_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "cache" / "reference.bin"
    destination.parent.mkdir(parents=True)
    partial = destination.with_name("reference.bin.part")
    partial.write_bytes(b"complete")
    destination.with_name("reference.bin.part.json").write_text(
        json.dumps(
            {
                "url": "https://example.test/reference.bin",
                "etag": '"reference-v1"',
                "last_modified": None,
                "total_bytes": 8,
            }
        )
    )
    reporter = ProvisioningReporter(
        tmp_path / "events.jsonl",
        stream=io.StringIO(),
        progress_interval_seconds=0,
    )

    def fake_urlopen(request: urllib.request.Request, *, timeout: int) -> _HTTPResponse:
        assert timeout == 120
        assert request.get_header("Range") == "bytes=8-"
        raise urllib.error.HTTPError(
            request.full_url,
            416,
            "range not satisfiable",
            {"Content-Range": "bytes */8"},
            None,
        )

    monkeypatch.setattr(progress_module.urllib.request, "urlopen", fake_urlopen)

    resumable_download(
        "https://example.test/reference.bin",
        destination,
        reporter,
        attempts=1,
        sleep=lambda _seconds: None,
    )

    assert destination.read_bytes() == b"complete"
    assert not partial.exists()


def test_resumable_download_rejects_nonzero_206_for_a_fresh_request(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "cache" / "reference.bin"
    reporter = ProvisioningReporter(
        tmp_path / "events.jsonl",
        stream=io.StringIO(),
        progress_interval_seconds=0,
    )
    responses = iter(
        (
            _HTTPResponse(
                b"wrong",
                status=206,
                headers={"Content-Length": "5", "Content-Range": "bytes 5-9/10"},
            ),
            _HTTPResponse(
                b"0123456789",
                status=200,
                headers={"Content-Length": "10"},
            ),
        )
    )
    requested_ranges: list[str | None] = []

    def fake_urlopen(request: urllib.request.Request, *, timeout: int) -> _HTTPResponse:
        assert timeout == 120
        requested_ranges.append(request.get_header("Range"))
        return next(responses)

    monkeypatch.setattr(progress_module.urllib.request, "urlopen", fake_urlopen)

    resumable_download(
        "https://example.test/reference.bin",
        destination,
        reporter,
        attempts=2,
        sleep=lambda _seconds: None,
    )

    assert requested_ranges == [None, None]
    assert destination.read_bytes() == b"0123456789"
