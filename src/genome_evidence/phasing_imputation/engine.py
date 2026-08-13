"""Narrow, safe Beagle 5.5 adapter. Analysis performs no network access."""

import json
import subprocess
from hashlib import sha256
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field


class BeagleEngine(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")
    jar: Path
    java: Path = Path("java")
    version: str = "5.5-27Feb25.75f"
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    byte_size: int = Field(gt=0)


def validate_beagle(engine: BeagleEngine) -> str:
    if not engine.jar.is_file() or engine.jar.stat().st_size != engine.byte_size:
        raise ValueError("Beagle JAR is missing or has the wrong byte size")
    if sha256(engine.jar.read_bytes()).hexdigest() != engine.sha256:
        raise ValueError("Beagle JAR checksum mismatch")
    result = subprocess.run(
        [str(engine.java), "-jar", str(engine.jar)], capture_output=True, text=True, timeout=30
    )
    identity = result.stdout + result.stderr
    if "beagle 5.5" not in identity.lower() or "27Feb25.75f" not in identity:
        raise ValueError("unexpected Java/Beagle identity")
    return identity.splitlines()[0] if identity.splitlines() else engine.version


def run_beagle(
    engine: BeagleEngine,
    *,
    target: Path,
    reference: Path,
    genetic_map: Path,
    output_prefix: Path,
    chromosome: str,
    threads: int,
    seed: int,
    ne: int,
    timeout: int,
    private_log: Path,
) -> Path:
    """Run validated local Beagle with explicit arguments and private captured logs."""
    validate_beagle(engine)
    for path in (target, reference, genetic_map, output_prefix):
        if any(ord(c) < 32 for c in str(path)):
            raise ValueError("control character in engine path")
    command = [
        str(engine.java),
        f"-Xmx{2048}m",
        "-jar",
        str(engine.jar),
        f"gt={target}",
        f"ref={reference}",
        f"map={genetic_map}",
        f"out={output_prefix}",
        f"chrom={chromosome}",
        f"nthreads={threads}",
        f"seed={seed}",
        f"ne={ne}",
    ]
    result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    private_log.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(
            json.dumps({"engine": "beagle", "status": "failed", "returncode": result.returncode})
        )
    output = output_prefix.with_suffix(".vcf.gz")
    if not output.is_file():
        raise RuntimeError("Beagle completed without its declared VCF output")
    return output
