import json
from pathlib import Path

from genome_evidence.workspace import (
    WorkspaceConfig,
    import_23andme_source,
    initialize_workspace,
    resolve_personal_m2_run,
    run_personal_m1_m2,
)


def test_personal_workflow_publishes_and_resolves_m2(tmp_path: Path) -> None:
    root = initialize_workspace(tmp_path / "private", WorkspaceConfig(subject_id="subject-0001"))
    inbox = root / "inputs/inbox/23andme/source.txt"
    inbox.write_text(
        "# genome build: GRCh38\nsynthetic_ref\t1\t5\tAA\nsynthetic_no_call\t1\t8\t--\n"
    )
    import_23andme_source(root, inbox, "subject-0001")
    marker_directory = root / "references/markers/23andme"
    (marker_directory / "marker-definitions.json").write_text(
        json.dumps(
            [
                {
                    "marker_id": "synthetic_ref",
                    "assembly": "GRCh38",
                    "chromosome": "1",
                    "position": 5,
                    "reference": "A",
                    "alternate": "G",
                    "orientation": "none",
                    "orientation_authoritative": True,
                },
                {
                    "marker_id": "synthetic_no_call",
                    "assembly": "GRCh38",
                    "chromosome": "1",
                    "position": 8,
                    "reference": "A",
                    "alternate": "T",
                    "orientation": "none",
                    "orientation_authoritative": True,
                },
            ]
        )
    )
    (root / "references/genome/grch38/GRCh38.fa").write_text(">1\n" + "A" * 20 + "\n")

    result = run_personal_m1_m2(
        root,
        "subject-0001",
        {},
        working_root=tmp_path / "ephemeral",
    )

    assert result.observation_count == 2
    assert result.mapping_count == 2
    assert result.canonical_genotype_count == 1
    assert result.m1_run.parent == root / "runs/m1_ingestion"
    assert result.m2_run.parent == root / "runs/m2_normalization"
    assert resolve_personal_m2_run(root, {}) == result.m2_run
    assert (
        resolve_personal_m2_run(
            root, {"GENOME_EVIDENCE_NORMALIZATION_RUN": str(result.m2_run.relative_to(root))}
        )
        == result.m2_run
    )
