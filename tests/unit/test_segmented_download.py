import io
import threading
from email.message import Message
from pathlib import Path
from urllib.request import Request

import pytest

from genome_evidence.workspace.provisioning_progress import ProvisioningReporter
from genome_evidence.workspace.segmented_download import RemoteIdentity, segmented_download


class _Response:
    def __init__(self, payload: bytes, start: int, end: int, total: int) -> None:
        self.status = 206
        self._stream = io.BytesIO(payload)
        self.headers = Message()
        self.headers["Content-Range"] = f"bytes {start}-{end}/{total}"
        self.headers["ETag"] = '"synthetic-v1"'

    def __enter__(self) -> "_Response":
        return self

    def __exit__(self, *_args: object) -> None:
        return None

    def read(self, size: int = -1) -> bytes:
        return self._stream.read(size)

    def getcode(self) -> int:
        return self.status


def test_segmented_download_reconstructs_with_bounded_overlap(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Explicitly synthetic bytes; no UCSC or private source is contacted.
    payload = bytes(range(251)) * 20
    lock = threading.Lock()
    barrier = threading.Barrier(2)
    active = 0
    maximum = 0

    def open_request(request: Request, **_kwargs: object) -> _Response:
        nonlocal active, maximum
        start, end = map(int, request.headers["Range"].removeprefix("bytes=").split("-"))
        with lock:
            active += 1
            maximum = max(maximum, active)
        if start < 2048:
            barrier.wait(timeout=2)
        response = _Response(payload[start : end + 1], start, end, len(payload))
        with lock:
            active -= 1
        return response

    monkeypatch.setattr(
        "genome_evidence.workspace.segmented_download.urllib.request.urlopen", open_request
    )
    reporter = ProvisioningReporter(tmp_path / "events.jsonl", stream=io.StringIO())
    destination, manifest = segmented_download(
        "https://example.invalid/synthetic.bb",
        tmp_path / "synthetic.bb",
        reporter,
        identity=RemoteIdentity(
            "https://example.invalid/synthetic.bb", len(payload), '"synthetic-v1"', None
        ),
        concurrency=3,
        segment_bytes=1024,
        sleep=lambda _seconds: None,
    )

    assert destination.read_bytes() == payload
    assert maximum >= 2
    assert manifest["byte_size"] == len(payload)
    assert destination.with_name("synthetic.bb.COMPLETED.json").is_file()
    assert "rs987654321" not in reporter.log_path.read_text(encoding="utf-8")
