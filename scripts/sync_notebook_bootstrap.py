"""Synchronize the unavoidable canonical pre-import cell across notebooks."""

import json
from pathlib import Path

BOOTSTRAP = """import importlib
import importlib.metadata
import json
import subprocess
from hashlib import sha256
from pathlib import Path

if PROFILE not in {"personal_drive", "synthetic_ci"}:
    raise ValueError("PROFILE must be personal_drive or synthetic_ci")

CHECKOUT = Path("/content/genome-evidence-src")
if PROFILE == "personal_drive":
    if "google.colab" in sys.modules:
        from google.colab import drive

        drive.mount("/content/drive", force_remount=False)
    if CHECKOUT.exists():
        remote = subprocess.run(
            ["git", "-C", str(CHECKOUT), "remote", "get-url", "origin"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout.strip()
        if remote != REPOSITORY_URL:
            raise RuntimeError("Unexpected checkout remote; move the checkout aside and rerun")
        dirty = subprocess.run(
            ["git", "-C", str(CHECKOUT), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout
        if dirty:
            raise RuntimeError("Checkout is dirty; preserve or move it aside and rerun")
    else:
        subprocess.run(
            ["git", "clone", "--no-checkout", REPOSITORY_URL, str(CHECKOUT)],
            check=True,
            timeout=180,
        )
    subprocess.run(
        ["git", "-C", str(CHECKOUT), "fetch", "--force", "origin", REPOSITORY_REF],
        check=True,
        timeout=180,
    )
    RESOLVED_COMMIT = subprocess.run(
        ["git", "-C", str(CHECKOUT), "rev-parse", "--verify", "FETCH_HEAD^{commit}"],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    ).stdout.strip()
    subprocess.run(
        ["git", "-C", str(CHECKOUT), "checkout", "--detach", RESOLVED_COMMIT],
        check=True,
        timeout=60,
    )
    previously_imported = sys.modules.get("genome_evidence")
    if previously_imported is not None:
        previous_file = Path(getattr(previously_imported, "__file__", "")).resolve()
        if not previous_file.is_relative_to(CHECKOUT.resolve()):
            raise RuntimeError(
                "genome_evidence was already imported elsewhere; restart the runtime"
            )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            "-e",
            f"{CHECKOUT}[notebook]",
        ],
        check=True,
        timeout=600,
    )
else:
    RESOLVED_COMMIT = "installed-ci-package"

genome_evidence = importlib.import_module("genome_evidence")
PACKAGE_ORIGIN = Path(genome_evidence.__file__).resolve()
if PROFILE == "personal_drive" and not PACKAGE_ORIGIN.is_relative_to(CHECKOUT.resolve()):
    raise RuntimeError("genome_evidence import origin is outside the resolved checkout")
INSTALLED_VERSION = importlib.metadata.version("genome-evidence")
LOCK_SHA256 = (
    sha256((CHECKOUT / "uv.lock").read_bytes()).hexdigest() if PROFILE == "personal_drive" else None
)
SANITIZED_IMPORT_PATH = (
    str(PACKAGE_ORIGIN.relative_to(CHECKOUT))
    if PROFILE == "personal_drive"
    else "installed-ci-package"
)
BOOTSTRAP_STATUS = {
    "profile": PROFILE,
    "requested_ref": REPOSITORY_REF,
    "resolved_commit": RESOLVED_COMMIT,
    "version": INSTALLED_VERSION,
    "import_path": SANITIZED_IMPORT_PATH,
    "lock_sha256": LOCK_SHA256,
    "lock_equivalent": PROFILE != "personal_drive",
}
print(json.dumps(BOOTSTRAP_STATUS, sort_keys=True))"""


def synchronize() -> None:
    for path in sorted(Path("notebooks").glob("*.ipynb")):
        notebook = json.loads(path.read_text())
        tagged = [
            cell
            for cell in notebook["cells"]
            if "genome-evidence-bootstrap" in cell.get("metadata", {}).get("tags", [])
        ]
        if len(tagged) != 1:
            raise ValueError(f"{path}: expected exactly one bootstrap cell")
        tagged[0]["source"] = BOOTSTRAP.splitlines(keepends=True)
        path.write_text(json.dumps(notebook, indent=1) + "\n")


if __name__ == "__main__":
    synchronize()
