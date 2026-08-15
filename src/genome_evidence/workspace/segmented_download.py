"""Identity-bound, resumable segmented HTTP downloads for large public resources."""

import http.client
import json
import os
import random
import re
import shutil
import threading
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import suppress
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path

from .provisioning_progress import MiB, ProvisioningReporter, resumable_download

_RANGE = re.compile(r"bytes (\d+)-(\d+)/(\d+)")


class RangeUnsupported(RuntimeError):
    """The origin ignored a byte range, requiring single-stream fallback."""


@dataclass(frozen=True)
class RemoteIdentity:
    canonical_url: str
    total_bytes: int
    etag: str | None
    last_modified: str | None

    @property
    def key(self) -> str:
        return sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()


def remote_identity(url: str, *, timeout_seconds: int = 120) -> RemoteIdentity:
    """Resolve canonical URL, exact size, and available HTTP validators."""
    request = urllib.request.Request(
        url, method="HEAD", headers={"User-Agent": "genome-evidence/0.4.2"}
    )
    try:
        response = urllib.request.urlopen(request, timeout=timeout_seconds)
    except urllib.error.HTTPError as error:
        if error.code not in {400, 403, 405, 501}:
            raise
        response = urllib.request.urlopen(
            urllib.request.Request(
                url,
                headers={"Range": "bytes=0-0", "User-Agent": "genome-evidence/0.4.2"},
            ),
            timeout=timeout_seconds,
        )
    with response:
        content_range = response.headers.get("Content-Range", "")
        match = _RANGE.fullmatch(content_range)
        length = int(match.group(3)) if match else int(response.headers.get("Content-Length", "0"))
        if length <= 0:
            raise OSError("remote object has no exact positive byte length")
        return RemoteIdentity(
            response.geturl(),
            length,
            response.headers.get("ETag"),
            response.headers.get("Last-Modified"),
        )


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".writing")
    temporary.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")
    with temporary.open("rb") as handle:
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _segment(
    identity: RemoteIdentity,
    start: int,
    end: int,
    directory: Path,
    reporter: ProvisioningReporter,
    cancel: threading.Event,
    *,
    attempts: int,
    timeout_seconds: int,
    chunk_bytes: int,
    sleep: Callable[[float], None],
    jitter: Callable[[], float],
) -> Path:
    key = f"{start:020d}-{end:020d}"
    partial, state_path, completed_path = (
        directory / f"{key}.part",
        directory / f"{key}.json",
        directory / f"{key}.COMPLETED.json",
    )
    expected = end - start + 1
    state = {}
    with suppress(OSError, UnicodeError, json.JSONDecodeError):
        state = json.loads(state_path.read_text(encoding="utf-8"))
    compatible = (
        state.get("identity_key") == identity.key
        and state.get("start") == start
        and state.get("end") == end
    )
    if not compatible:
        partial.unlink(missing_ok=True)
        completed_path.unlink(missing_ok=True)
    if completed_path.is_file() and partial.is_file() and compatible:
        try:
            completed = json.loads(completed_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            completed = {}
        if partial.stat().st_size == expected and completed.get("sha256") == _hash(partial):
            return partial
    _write_json(
        state_path,
        {"identity_key": identity.key, "start": start, "end": end, "expected_length": expected},
    )
    last: BaseException | None = None
    for attempt in range(1, attempts + 1):
        if cancel.is_set():
            raise RangeUnsupported("another segment observed an origin without range support")
        offset = partial.stat().st_size if partial.is_file() else 0
        if offset > expected:
            partial.unlink()
            offset = 0
        requested_start = start + offset
        headers = {"Range": f"bytes={requested_start}-{end}", "User-Agent": "genome-evidence/0.4.2"}
        validator = identity.etag or identity.last_modified
        if validator:
            headers["If-Range"] = validator
        try:
            with urllib.request.urlopen(
                urllib.request.Request(identity.canonical_url, headers=headers),
                timeout=timeout_seconds,
            ) as response:
                status = getattr(response, "status", response.getcode())
                if status == 200:
                    cancel.set()
                    raise RangeUnsupported("origin ignored HTTP Range")
                match = _RANGE.fullmatch(response.headers.get("Content-Range", ""))
                if (
                    status != 206
                    or match is None
                    or tuple(map(int, match.groups()))
                    != (
                        requested_start,
                        end,
                        identity.total_bytes,
                    )
                ):
                    raise OSError("origin returned an incompatible Content-Range")
                if identity.etag and response.headers.get("ETag") not in {None, identity.etag}:
                    raise OSError("remote ETag changed during segmented transfer")
                if identity.last_modified and response.headers.get("Last-Modified") not in {
                    None,
                    identity.last_modified,
                }:
                    raise OSError("remote Last-Modified changed during segmented transfer")
                with partial.open("ab" if offset else "wb") as output:
                    while data := response.read(chunk_bytes):
                        output.write(data)
                    output.flush()
                    os.fsync(output.fileno())
            if partial.stat().st_size != expected:
                raise OSError("ranged response ended before the requested segment was complete")
            digest = _hash(partial)
            _write_json(
                completed_path,
                {
                    "identity_key": identity.key,
                    "start": start,
                    "end": end,
                    "byte_size": expected,
                    "sha256": digest,
                },
            )
            return partial
        except urllib.error.HTTPError as error:
            if error.code == 416 and partial.is_file() and partial.stat().st_size == expected:
                digest = _hash(partial)
                _write_json(
                    completed_path,
                    {
                        "identity_key": identity.key,
                        "start": start,
                        "end": end,
                        "byte_size": expected,
                        "sha256": digest,
                    },
                )
                return partial
            last = error
        except RangeUnsupported:
            raise
        except (OSError, TimeoutError, urllib.error.URLError, http.client.HTTPException) as error:
            last = error
        reporter.warning(
            "common.download.retry",
            "A common-file segment was interrupted; retained bytes will resume.",
            attempt=attempt,
            attempts=attempts,
            segment_start=start,
            segment_end=end,
            error_type=type(last).__name__,
        )
        if attempt < attempts:
            sleep(min(2 ** (attempt - 1), 30) + jitter())
    raise RuntimeError("segmented transfer exhausted bounded retries") from last


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(MiB), b""):
            digest.update(chunk)
    return digest.hexdigest()


def segmented_download(
    url: str,
    destination: Path,
    reporter: ProvisioningReporter,
    *,
    concurrency: int = 8,
    segment_bytes: int = 64 * MiB,
    attempts: int = 5,
    timeout_seconds: int = 120,
    chunk_bytes: int = 4 * MiB,
    sleep: Callable[[float], None] = time.sleep,
    jitter: Callable[[], float] = random.random,
    validate: Callable[[Path], None] | None = None,
    identity: RemoteIdentity | None = None,
) -> tuple[Path, dict[str, object]]:
    """Download ranges concurrently, assemble locally, validate, then commit a manifest last."""
    if not 1 <= concurrency <= 12:
        raise ValueError("segment concurrency must be between 1 and 12")
    identity = identity or remote_identity(url, timeout_seconds=timeout_seconds)
    completion_path = destination.with_name(destination.name + ".COMPLETED.json")
    if destination.is_file() and completion_path.is_file():
        try:
            installed = json.loads(completion_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            installed = {}
        if (
            installed.get("schema") == "genome-evidence-segmented-download/v1"
            and installed.get("remote_identity") == asdict(identity)
            and installed.get("byte_size") == identity.total_bytes
            and installed.get("sha256") == _hash(destination)
        ):
            reporter.success(
                "common.download.reuse",
                "Reused a completed destination after explicit manifest verification; "
                "no segment assembly or download was performed.",
                reused_complete_segments=installed.get("segment_count", 0),
                session_downloaded_bytes=0,
                final_file_verification=True,
                byte_size=identity.total_bytes,
            )
            return destination, installed
    directory = destination.with_name(destination.name + ".segments")
    directory.mkdir(parents=True, exist_ok=True)
    ranges = [
        (start, min(start + segment_bytes - 1, identity.total_bytes - 1))
        for start in range(0, identity.total_bytes, segment_bytes)
    ]
    cancel = threading.Event()
    started = time.monotonic()
    retained_bytes = 0
    reused_segments = 0
    for start, end in ranges:
        key = f"{start:020d}-{end:020d}"
        partial = directory / f"{key}.part"
        state_path = directory / f"{key}.json"
        completed_path = directory / f"{key}.COMPLETED.json"
        try:
            state = json.loads(state_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            state = {}
        if state.get("identity_key") == identity.key and partial.is_file():
            retained_bytes += min(partial.stat().st_size, end - start + 1)
            if completed_path.is_file() and partial.stat().st_size == end - start + 1:
                reused_segments += 1
    try:
        with ThreadPoolExecutor(max_workers=concurrency, thread_name_prefix="common-range") as pool:
            futures = {
                pool.submit(
                    _segment,
                    identity,
                    start,
                    end,
                    directory,
                    reporter,
                    cancel,
                    attempts=attempts,
                    timeout_seconds=timeout_seconds,
                    chunk_bytes=chunk_bytes,
                    sleep=sleep,
                    jitter=jitter,
                ): index
                for index, (start, end) in enumerate(ranges)
            }
            outputs: list[Path | None] = [None] * len(ranges)
            completed = 0
            for future in as_completed(futures):
                outputs[futures[future]] = future.result()
                completed += outputs[futures[future]].stat().st_size  # type: ignore[union-attr]
                reporter.progress(
                    "common.download.progress",
                    destination.name,
                    completed,
                    identity.total_bytes,
                    unit="bytes",
                    started_at=started,
                    initial_completed=retained_bytes,
                    force=True,
                    completed_segments=sum(item is not None for item in outputs),
                    total_segments=len(outputs),
                    configured_workers=concurrency,
                    session_downloaded_bytes=max(completed - retained_bytes, 0),
                    resumed_partial_bytes=retained_bytes,
                    reused_complete_segments=reused_segments,
                )
    except RangeUnsupported:
        cancel.set()
        reporter.warning(
            "common.download.single-stream",
            "The origin ignored Range; switching safely to one resumable stream.",
            url=identity.canonical_url,
        )
        result = resumable_download(
            identity.canonical_url,
            destination,
            reporter,
            attempts=attempts,
            sleep=sleep,
            timeout_seconds=timeout_seconds,
            chunk_bytes=chunk_bytes,
        )
        if result.stat().st_size != identity.total_bytes:
            raise OSError(
                "single-stream fallback length differs from probed remote identity"
            ) from None
    else:
        reporter.info(
            "common.download.assemble",
            "ASSEMBLE verified segments into the destination file.",
            segment_count=len(ranges),
        )
        temporary = destination.with_name(destination.name + ".assembling")
        with temporary.open("wb") as target:
            for output in outputs:
                if output is None:
                    raise RuntimeError("segment output is missing")
                with output.open("rb") as source:
                    shutil.copyfileobj(source, target, length=4 * MiB)
            target.flush()
            os.fsync(target.fileno())
        if temporary.stat().st_size != identity.total_bytes:
            raise OSError("assembled common file has the wrong byte length")
        os.replace(temporary, destination)
    if validate:
        reporter.info("common.download.verify", "VERIFY the assembled destination file.")
        validate(destination)
    manifest = {
        "schema": "genome-evidence-segmented-download/v1",
        "remote_identity": asdict(identity),
        "sha256": _hash(destination),
        "byte_size": destination.stat().st_size,
        "segment_count": len(ranges),
        "segment_concurrency": concurrency,
    }
    _write_json(completion_path, manifest)
    reporter.success(
        "common.download.complete",
        "Completed and verified a common dbSNP file.",
        byte_size=identity.total_bytes,
        segment_count=len(ranges),
        sha256=manifest["sha256"],
    )
    return destination, manifest
