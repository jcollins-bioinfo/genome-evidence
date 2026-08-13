"""Shared, injectable notebook profile and exact-source bootstrap support."""

import importlib
import importlib.metadata
import os
import subprocess
import sys
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path

PROFILES = ("personal_drive", "synthetic_ci")
REPOSITORY_URL = "https://github.com/jcollins-bioinfo/genome-evidence.git"
DRIVE_MOUNT = Path("/content/drive")
DEFAULT_WORKSPACE = Path("/content/drive/MyDrive/genome-evidence-private")
DEFAULT_CHECKOUT = Path("/content/genome-evidence-src")


@dataclass(frozen=True)
class NotebookSettings:
    profile: str
    repository_url: str
    repository_ref: str
    workspace_root: Path
    subject_id: str


def resolve_settings(
    environment: Mapping[str, str] | None = None, *, in_colab: bool | None = None
) -> NotebookSettings:
    """Resolve exactly one profile; personal mode never silently becomes synthetic."""
    env = os.environ if environment is None else environment
    colab = "google.colab" in sys.modules if in_colab is None else in_colab
    profile = env.get("GENOME_EVIDENCE_PROFILE", "personal_drive" if colab else "synthetic_ci")
    if profile not in PROFILES:
        raise ValueError(f"GENOME_EVIDENCE_PROFILE must be one of {PROFILES}")
    return NotebookSettings(
        profile,
        REPOSITORY_URL,
        env.get("GENOME_EVIDENCE_GIT_REF", "main"),
        Path(env.get("GENOME_EVIDENCE_WORKSPACE", str(DEFAULT_WORKSPACE))),
        env.get("GENOME_EVIDENCE_SUBJECT_ID", "subject-0001"),
    )


def prepare_personal_runtime(
    settings: NotebookSettings,
    *,
    checkout: Path = DEFAULT_CHECKOUT,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
) -> dict[str, str]:
    """Safely reuse/clone, detach, install, and prove the requested public checkout."""
    if settings.profile != "personal_drive":
        raise ValueError("source installation is forbidden outside personal_drive")
    if checkout.exists():
        remote = runner(
            ["git", "-C", str(checkout), "remote", "get-url", "origin"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout.strip()
        if remote != settings.repository_url:
            raise ValueError("wrong source remote; move the checkout aside and rerun")
        dirty = runner(
            ["git", "-C", str(checkout), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout
        if dirty:
            raise ValueError("source checkout is dirty; preserve or move it aside and rerun")
        runner(["git", "-C", str(checkout), "fetch", "--tags", "origin"], check=True, timeout=120)
    else:
        runner(
            ["git", "clone", "--no-checkout", settings.repository_url, str(checkout)],
            check=True,
            timeout=180,
        )
    resolved = runner(
        [
            "git",
            "-C",
            str(checkout),
            "rev-parse",
            "--verify",
            f"{settings.repository_ref}^{{commit}}",
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    ).stdout.strip()
    runner(["git", "-C", str(checkout), "checkout", "--detach", resolved], check=True, timeout=60)
    runner(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            "-e",
            f"{checkout}[notebook]",
        ],
        check=True,
        timeout=600,
    )
    existing = sys.modules.get("genome_evidence")
    if (
        existing is not None
        and existing.__file__ is not None
        and not Path(existing.__file__).resolve().is_relative_to(checkout.resolve())
    ):
        raise RuntimeError("genome_evidence was imported from another location; restart runtime")
    package = importlib.import_module("genome_evidence")
    if package.__file__ is None:
        raise RuntimeError("installed package has no inspectable file origin")
    origin = Path(package.__file__).resolve()
    if not origin.is_relative_to(checkout.resolve()):
        raise RuntimeError("installed package origin is outside the resolved checkout")
    return {
        "commit": resolved,
        "version": importlib.metadata.version("genome-evidence"),
        "import_path": str(origin.relative_to(checkout.resolve())),
    }
