"""Privacy-safe progress telemetry and resumable HTTP transfers for workspace provisioning."""

import http.client
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
import uuid
from collections.abc import Callable
from contextlib import suppress
from datetime import UTC, datetime
from pathlib import Path
from types import TracebackType
from typing import IO, Any, Self

MiB = 1024 * 1024
Clock = Callable[[], float]
Sleeper = Callable[[float], None]


def _format_amount(value: float, unit: str) -> str:
    if unit == "bytes":
        if value >= 1024**3:
            return f"{value / 1024**3:.2f} GiB"
        if value >= MiB:
            return f"{value / MiB:.1f} MiB"
        if value >= 1024:
            return f"{value / 1024:.1f} KiB"
        return f"{value:.0f} B"
    return f"{value:,.0f} {unit}"


def _format_duration(seconds: float | None) -> str:
    if seconds is None or seconds < 0:
        return "--:--"
    rounded = int(seconds + 0.5)
    hours, remainder = divmod(rounded, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}" if hours else f"{minutes:02d}:{seconds:02d}"


class ProvisioningReporter:
    """Append structured JSONL events while mirroring concise progress to Jupyter stdout."""

    def __init__(
        self,
        log_path: Path,
        *,
        stream: IO[str] | None = None,
        clock: Clock = time.monotonic,
        progress_interval_seconds: float = 5.0,
    ) -> None:
        self.log_path = log_path
        self.session_id = uuid.uuid4().hex[:12]
        self.stream = sys.stdout if stream is None else stream
        self._clock = clock
        self._interval = progress_interval_seconds
        self._last_progress: dict[str, float] = {}
        self._log_available = True
        self._log_error_reported = False
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.touch(exist_ok=True)
        except OSError:
            self._log_available = False
        else:
            with suppress(OSError):
                os.chmod(log_path, 0o600)

    def event(self, level: str, event: str, message: str, **fields: Any) -> None:
        timestamp = datetime.now(UTC).isoformat()
        record = {
            "timestamp": timestamp,
            "session_id": self.session_id,
            "level": level,
            "event": event,
            "message": message,
            **fields,
        }
        clock = timestamp[11:19]
        if not self._log_available:
            try:
                self.log_path.parent.mkdir(parents=True, exist_ok=True)
                self.log_path.touch(exist_ok=True)
                self._log_available = True
            except OSError:
                pass
        if self._log_available:
            try:
                with self.log_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(record, sort_keys=True, default=str) + "\n")
                    handle.flush()
            except OSError:
                self._log_available = False
        if not self._log_available and not self._log_error_reported:
            self._log_error_reported = True
            print(
                f"[00B {clock}] WARN  log.unavailable          "
                "The durable log is temporarily unavailable; stdout progress will continue.",
                file=self.stream,
                flush=True,
            )
        elif self._log_available and self._log_error_reported:
            self._log_error_reported = False
            print(
                f"[00B {clock}] DONE  log.recovered            Durable JSONL logging resumed.",
                file=self.stream,
                flush=True,
            )
        print(f"[00B {clock}] {level:<5} {event:<24} {message}", file=self.stream, flush=True)

    def info(self, event: str, message: str, **fields: Any) -> None:
        self.event("INFO", event, message, **fields)

    def success(self, event: str, message: str, **fields: Any) -> None:
        self.event("DONE", event, message, **fields)

    def warning(self, event: str, message: str, **fields: Any) -> None:
        self.event("WARN", event, message, **fields)

    def error(self, event: str, message: str, **fields: Any) -> None:
        self.event("ERROR", event, message, **fields)

    def progress(
        self,
        event: str,
        label: str,
        completed: int,
        total: int | None,
        *,
        unit: str,
        started_at: float,
        initial_completed: int = 0,
        force: bool = False,
        **fields: Any,
    ) -> None:
        now = self._clock()
        if not force and now - self._last_progress.get(event, float("-inf")) < self._interval:
            return
        self._last_progress[event] = now
        elapsed = max(now - started_at, 1e-9)
        session_completed = max(completed - initial_completed, 0)
        rate = session_completed / elapsed
        eta = ((total - completed) / rate) if total is not None and rate > 0 else None
        proportion = completed / total if total else None
        total_text = _format_amount(total, unit) if total is not None else "unknown"
        rate_text = f"{_format_amount(rate, unit)}/s"
        percent_text = f"{proportion * 100:5.1f}%" if proportion is not None else "  n/a"
        message = (
            f"{label} | {_format_amount(completed, unit)} / {total_text} | "
            f"{percent_text} | {rate_text} | ETA {_format_duration(eta)}"
        )
        self.event(
            "PROG",
            event,
            message,
            completed=completed,
            total=total,
            unit=unit,
            rate_per_second=rate,
            session_completed=session_completed,
            elapsed_seconds=elapsed,
            eta_seconds=eta,
            **fields,
        )

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        if exception is not None:
            self.error(
                "session.paused",
                "Provisioning paused; durable checkpoints were retained. "
                "Rerun this cell to resume.",
                exception_type=type(exception).__name__,
            )


def resumable_download(
    url: str,
    destination: Path,
    reporter: ProvisioningReporter,
    *,
    attempts: int = 5,
    sleep: Sleeper = time.sleep,
    timeout_seconds: int = 120,
    chunk_bytes: int = 4 * MiB,
) -> Path:
    """Download with HTTP Range continuation, bounded retries, rates, and durable partial bytes."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file():
        reporter.info(
            "download.reuse",
            f"Reusing completed cache file {destination.name} "
            f"({_format_amount(destination.stat().st_size, 'bytes')})",
            url=url,
            path=str(destination),
            byte_size=destination.stat().st_size,
        )
        return destination
    partial = destination.with_name(destination.name + ".part")
    state_path = destination.with_name(destination.name + ".part.json")
    last_error: BaseException | None = None
    for attempt in range(1, attempts + 1):
        offset = partial.stat().st_size if partial.is_file() else 0
        headers = {"User-Agent": "genome-evidence/0.4.1"}
        try:
            candidate_state = json.loads(state_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            candidate_state = {}
        state = candidate_state if isinstance(candidate_state, dict) else {}
        validator: str | None = None
        stored_total = state.get("total_bytes")
        stored_validator = state.get("etag") or state.get("last_modified")
        if offset and not (
            state.get("url") == url
            and isinstance(stored_total, int)
            and stored_total >= offset
            and isinstance(stored_validator, str)
            and stored_validator
        ):
            reporter.warning(
                "download.partial.reset",
                "The retained partial has no trustworthy remote identity; restarting "
                "this file from byte zero.",
                url=url,
                retained_bytes=offset,
            )
            partial.unlink(missing_ok=True)
            state_path.unlink(missing_ok=True)
            offset = 0
            state = {}
        if offset:
            headers["Range"] = f"bytes={offset}-"
            validator = str(stored_validator)
            headers["If-Range"] = validator
        reporter.info(
            "download.attempt",
            f"Starting {destination.name}; attempt {attempt}/{attempts}, "
            f"resume offset {_format_amount(offset, 'bytes')}",
            url=url,
            path=str(destination),
            attempt=attempt,
            resume_offset=offset,
        )
        started_at = time.monotonic()
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                status = getattr(response, "status", response.getcode())
                content_range = response.headers.get("Content-Range")
                content_length = response.headers.get("Content-Length")
                response_etag = response.headers.get("ETag")
                response_modified = response.headers.get("Last-Modified")
                response_validator = response_etag or response_modified
                range_match = (
                    re.fullmatch(r"bytes (\d+)-(\d+)/(\d+|\*)", content_range)
                    if content_range
                    else None
                )
                if status == 206 and (range_match is None or not range_match.group(3).isdigit()):
                    raise OSError("server returned partial content without a verifiable total size")
                if not offset and (
                    status not in {200, 206}
                    or (status == 206 and (range_match is None or range_match.group(1) != "0"))
                ):
                    partial.unlink(missing_ok=True)
                    state_path.unlink(missing_ok=True)
                    raise OSError("server returned an incompatible response for a fresh download")
                if offset and status == 200:
                    reporter.warning(
                        "download.restart",
                        "Server did not honor the Range request; restarting this file "
                        "from byte zero.",
                        url=url,
                        http_status=status,
                    )
                    offset = 0
                elif offset and (
                    status != 206
                    or range_match is None
                    or int(range_match.group(1)) != offset
                    or int(range_match.group(3)) != stored_total
                    or (validator is not None and response_validator not in {None, validator})
                ):
                    partial.unlink(missing_ok=True)
                    state_path.unlink(missing_ok=True)
                    raise OSError(
                        "server returned an incompatible partial-content response; "
                        "the local partial was reset"
                    )
                total: int | None = None
                if range_match is not None:
                    total_text = range_match.group(3)
                    total = int(total_text) if total_text.isdigit() else None
                elif content_length and content_length.isdigit():
                    total = int(content_length) + (offset if status == 206 else 0)
                state_payload = {
                    "url": url,
                    "etag": response_etag,
                    "last_modified": response_modified,
                    "total_bytes": total,
                }
                state_temporary = state_path.with_name(state_path.name + ".tmp")
                state_temporary.write_text(
                    json.dumps(state_payload, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
                os.replace(state_temporary, state_path)
                mode = "ab" if offset and status == 206 else "wb"
                transferred = offset
                with partial.open(mode) as output:
                    while chunk := response.read(chunk_bytes):
                        output.write(chunk)
                        transferred += len(chunk)
                        reporter.progress(
                            f"download.{destination.name}",
                            destination.name,
                            transferred,
                            total,
                            unit="bytes",
                            started_at=started_at,
                            initial_completed=offset,
                            url=url,
                            attempt=attempt,
                        )
                    output.flush()
                    os.fsync(output.fileno())
                if total is not None and transferred != total:
                    raise OSError(
                        f"HTTP response ended at byte {transferred:,}; expected {total:,} bytes"
                    )
                os.replace(partial, destination)
                state_path.unlink(missing_ok=True)
                reporter.progress(
                    f"download.{destination.name}",
                    destination.name,
                    destination.stat().st_size,
                    total or destination.stat().st_size,
                    unit="bytes",
                    started_at=started_at,
                    initial_completed=offset,
                    force=True,
                    url=url,
                    attempt=attempt,
                )
                reporter.success(
                    "download.complete",
                    f"Completed {destination.name} "
                    f"({_format_amount(destination.stat().st_size, 'bytes')})",
                    url=url,
                    path=str(destination),
                    byte_size=destination.stat().st_size,
                )
                return destination
        except (
            OSError,
            TimeoutError,
            urllib.error.URLError,
            http.client.HTTPException,
        ) as error:
            if isinstance(error, urllib.error.HTTPError) and error.code == 416:
                content_range = error.headers.get("Content-Range") if error.headers else None
                complete_match = (
                    re.fullmatch(r"bytes \*/(\d+)", content_range) if content_range else None
                )
                expected = int(complete_match.group(1)) if complete_match else None
                retained = partial.stat().st_size if partial.is_file() else 0
                if expected is not None and retained == expected:
                    os.replace(partial, destination)
                    state_path.unlink(missing_ok=True)
                    reporter.success(
                        "download.complete",
                        f"Promoted the complete resumed file {destination.name} "
                        f"({_format_amount(retained, 'bytes')}).",
                        url=url,
                        path=str(destination),
                        byte_size=retained,
                        http_status=416,
                    )
                    return destination
                if expected is not None:
                    partial.unlink(missing_ok=True)
                    state_path.unlink(missing_ok=True)
                    reporter.warning(
                        "download.restart",
                        "The retained byte count does not match the remote object; "
                        "restarting this file from byte zero.",
                        url=url,
                        retained_bytes=retained,
                        remote_bytes=expected,
                        http_status=416,
                    )
            last_error = error
            retained = partial.stat().st_size if partial.is_file() else 0
            reporter.warning(
                "download.retry",
                f"Transfer interrupted; retained {_format_amount(retained, 'bytes')} for resume.",
                url=url,
                attempt=attempt,
                attempts=attempts,
                retained_bytes=retained,
                error_type=type(error).__name__,
                error_message=str(error)[:1000],
            )
            if attempt < attempts:
                delay = min(2 ** (attempt - 1), 30)
                reporter.info(
                    "download.backoff",
                    f"Retrying {destination.name} in {delay} seconds.",
                    attempt=attempt,
                    next_attempt=attempt + 1,
                    delay_seconds=delay,
                )
                sleep(delay)
    raise RuntimeError(
        f"download failed after {attempts} resumable attempts: {url}"
    ) from last_error
