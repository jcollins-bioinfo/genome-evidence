"""Synchronize the unavoidable canonical pre-import cell across notebooks."""

import json
from pathlib import Path

BOOTSTRAP = """from genome_evidence.notebook_support import resolve_settings
SETTINGS = resolve_settings()
PROFILE = SETTINGS.profile
REPOSITORY_URL = SETTINGS.repository_url
REPOSITORY_REF = SETTINGS.repository_ref
WORKSPACE_ROOT = SETTINGS.workspace_root
SUBJECT_ID = SETTINGS.subject_id
print({"profile": PROFILE, "requested_ref": REPOSITORY_REF})
"""


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
