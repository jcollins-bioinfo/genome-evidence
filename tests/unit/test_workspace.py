import json
from hashlib import sha256
from pathlib import Path

import pytest

from genome_evidence.workspace import (
    WorkspaceConfig,
    import_23andme_source,
    initialize_workspace,
    publish_completed_run,
    resolve_latest_compatible_run,
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


def test_completed_run_publication_and_compatible_resolution(tmp_path: Path) -> None:
    root = initialize_workspace(tmp_path / "private", WorkspaceConfig(subject_id="subject-0001"))
    source = tmp_path / "computed-m2"
    source.mkdir()
    artifact = source / "aggregate.json"
    artifact.write_text('{"mapped": 2}\n')
    run_id = "synthetic-m2-run"
    (source / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "run_id": run_id,
                "artifacts": {"aggregate.json": sha256(artifact.read_bytes()).hexdigest()},
            }
        )
    )

    published = publish_completed_run(
        root,
        source,
        "M2",
        {"target_assembly": "GRCh38", "source_sha256": "a" * 64},
    )

    assert published == root / "runs/m2_normalization" / run_id
    assert (published / "COMPLETED.json").is_file()
    assert not (root / "runs/m2_normalization/_incomplete" / run_id).exists()
    assert resolve_latest_compatible_run(root, "M2", {"target_assembly": "GRCh38"}) == published
    with pytest.raises(ValueError, match="incompatible"):
        resolve_latest_compatible_run(root, "M2", {"target_assembly": "GRCh37"})

    (published / "aggregate.json").write_text("tampered\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        resolve_latest_compatible_run(root, "M2", {"target_assembly": "GRCh38"})


def test_completed_run_rejects_undeclared_files(tmp_path: Path) -> None:
    root = initialize_workspace(tmp_path / "private", WorkspaceConfig(subject_id="subject-0001"))
    source = tmp_path / "computed-m1"
    source.mkdir()
    artifact = source / "observations.parquet"
    artifact.write_bytes(b"synthetic-only")
    (source / "manifest.json").write_text(
        json.dumps(
            {
                "run_id": "synthetic-m1-run",
                "artifacts": {"observations.parquet": sha256(artifact.read_bytes()).hexdigest()},
            }
        )
    )
    (source / "undeclared.txt").write_text("must not be copied")

    with pytest.raises(ValueError, match="undeclared"):
        publish_completed_run(root, source, "M1", {"resolved_build": "GRCh38"})


def test_completed_run_accepts_m5_sized_artifact_identities(tmp_path: Path) -> None:
    root = initialize_workspace(tmp_path / "private", WorkspaceConfig(subject_id="subject-0001"))
    source = tmp_path / "computed-m5"
    source.mkdir()
    artifact = source / "population_structure_qc.json"
    artifact.write_text('{"projection_status": "not_projected"}\n')
    (source / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "run_id": "synthetic-m5-run",
                "artifacts": {
                    artifact.name: {
                        "sha256": sha256(artifact.read_bytes()).hexdigest(),
                        "byte_size": artifact.stat().st_size,
                    }
                },
            }
        )
    )

    published = publish_completed_run(
        root,
        source,
        "M5",
        {"target_assembly": "GRCh38", "reference_model_id": "synthetic-model"},
    )

    assert resolve_latest_compatible_run(root, "M5", {"target_assembly": "GRCh38"}) == published
