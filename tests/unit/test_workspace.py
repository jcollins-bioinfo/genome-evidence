import json
from pathlib import Path

import pytest

from genome_evidence.workspace import (
    WorkspaceConfig,
    import_23andme_source,
    initialize_workspace,
    validate_workspace,
)


def test_workspace_idempotent_and_content_addressed(tmp_path: Path) -> None:
    root = tmp_path / "private"
    config = WorkspaceConfig(subject_id="subject-0001")
    assert initialize_workspace(root, config) == root
    assert initialize_workspace(root, config) == root
    source = root / "inputs/inbox/23andme/upload.txt"
    source.write_bytes(b"# synthetic only\nrs1\t1\t1\tAA\n")
    imported = import_23andme_source(root, source, "subject-0001")
    assert len(imported.name) == 64
    assert import_23andme_source(root, source, "subject-0001") == imported
    manifest = json.loads((imported / "source_manifest.json").read_text())
    assert manifest["data_file"] == "genome.txt"
    assert "upload" not in json.dumps(manifest)


def test_workspace_fails_closed(tmp_path: Path) -> None:
    root = tmp_path / "private"
    initialize_workspace(root, WorkspaceConfig(subject_id="subject-0001"))
    (root / "inputs/inbox/23andme/a").write_text("a")
    (root / "inputs/inbox/23andme/b").write_text("b")
    with pytest.raises(ValueError, match="exactly one"):
        import_23andme_source(root, None, "subject-0001")
    (root / "config/workspace.json").write_text('{"schema":"future"}')
    with pytest.raises(ValueError, match="unknown"):
        validate_workspace(root)
