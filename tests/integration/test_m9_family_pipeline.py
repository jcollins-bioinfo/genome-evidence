"""Synthetic end-to-end M9 input integrity, privacy, and publication tests."""

import json
from hashlib import sha256
from pathlib import Path

import polars as pl
from typer.testing import CliRunner

from genome_evidence.cli import app
from genome_evidence.family_analysis import analyze_family
from genome_evidence.workspace import WorkspaceConfig, initialize_workspace, publish_completed_run


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _m2(root: Path, number: int, calls: list[tuple[str, tuple[str, ...]]]) -> Path:
    workspace = initialize_workspace(root, WorkspaceConfig(subject_id=f"subject-{number:04d}"))
    run = root.parent / f"fabricated-staging-{number}"
    run.mkdir()
    variants = [
        {
            "variant_id": variant,
            "assembly": "GRCh38",
            "chromosome": "1",
            "position": index + 100,
            "reference": "A",
            "alternate": "G",
            "rsid": None,
        }
        for index, variant in enumerate(("fabricated-v1", "fabricated-v2", "fabricated-v3"))
    ]
    genotypes = [
        {
            "genotype_id": f"fabricated-gt-{number}-{variant}",
            "observation_reference": f"fabricated-observation-{number}-{variant}",
            "normalization_run_id": f"fabricated-m2-{number}",
            "variant_id": variant,
            "alleles": list(alleles),
            "ploidy": len(alleles),
            "call_status": "called",
        }
        for variant, alleles in calls
    ]
    pl.DataFrame(variants).write_parquet(run / "variants.parquet")
    pl.DataFrame(genotypes).write_parquet(run / "canonical_genotypes.parquet")
    pl.DataFrame().write_parquet(run / "observation_mappings.parquet")
    metadata = {
        "run_id": f"fabricated-m2-{number}",
        "m1_run_id": f"fabricated-m1-{number}",
        "target_assembly": "GRCh38",
        "algorithm": "m2-snv-1",
        "package_version": "fabricated",
        "resources": [],
    }
    (run / "normalization_metadata.json").write_text(json.dumps(metadata))
    artifacts = {path.name: _hash(path) for path in run.iterdir()}
    (run / "manifest.json").write_text(
        json.dumps(
            {"schema_version": 1, "run_id": f"fabricated-m2-{number}", "artifacts": artifacts}
        )
    )
    return publish_completed_run(workspace, run, "M2", {"target_assembly": "GRCh38"})


def test_family_pipeline_is_checksum_bound_atomic_and_cli_private(tmp_path: Path) -> None:
    parent_1 = _m2(
        tmp_path / "private-1",
        1,
        [
            ("fabricated-v1", ("A", "A")),
            ("fabricated-v2", ("A", "G")),
            ("fabricated-v3", ("A", "A")),
        ],
    )
    parent_2 = _m2(
        tmp_path / "private-2",
        2,
        [
            ("fabricated-v1", ("G", "G")),
            ("fabricated-v2", ("A", "G")),
            ("fabricated-v3", ("A", "A")),
        ],
    )
    child = _m2(
        tmp_path / "private-3",
        3,
        [
            ("fabricated-v1", ("A", "G")),
            ("fabricated-v2", ("A", "G")),
            ("fabricated-v3", ("G", "G")),
        ],
    )
    descriptor = tmp_path / "fabricated-pedigree.json"
    descriptor.write_text(
        json.dumps(
            {
                "schema": "genome-evidence-pedigree/v1",
                "family_id": "fabricated-family",
                "members": [
                    {
                        "member_id": "fabricated-parent-a",
                        "subject_id": "subject-0001",
                        "m2_run": str(parent_1),
                    },
                    {
                        "member_id": "fabricated-parent-b",
                        "subject_id": "subject-0002",
                        "m2_run": str(parent_2),
                    },
                    {
                        "member_id": "fabricated-child",
                        "subject_id": "subject-0003",
                        "m2_run": str(child),
                    },
                ],
                "relationships": [
                    {
                        "assertion_id": "fabricated-edge-a",
                        "parent_member_id": "fabricated-parent-a",
                        "child_member_id": "fabricated-child",
                        "source": "user_declared",
                    },
                    {
                        "assertion_id": "fabricated-edge-b",
                        "parent_member_id": "fabricated-parent-b",
                        "child_member_id": "fabricated-child",
                        "source": "user_declared",
                    },
                ],
            }
        )
    )
    output = tmp_path / "private-m9"
    result = analyze_family(descriptor, output)
    assert result.compatibility_counts == {"consistent": 2, "inconsistent": 1}
    assert (output / "COMPLETED.json").is_file()
    manifest = json.loads((output / "manifest.json").read_text())
    assert all(
        _hash(output / name) == identity["sha256"]
        for name, identity in manifest["artifacts"].items()
    )
    assert output.stat().st_mode & 0o077 == 0

    cli_descriptor = tmp_path / "fabricated-pedigree-cli.json"
    cli_descriptor.write_bytes(descriptor.read_bytes())
    cli = CliRunner().invoke(
        app, ["family", "validate-pedigree", "--pedigree", str(cli_descriptor)]
    )
    assert cli.exit_code == 0
    assert "fabricated-parent" not in cli.stdout and "subject-" not in cli.stdout


def test_tampered_m2_is_rejected(tmp_path: Path) -> None:
    run = _m2(tmp_path / "private-1", 1, [("fabricated-v1", ("A", "A"))])
    child_run = _m2(tmp_path / "private-2", 2, [("fabricated-v1", ("A", "A"))])
    (run / "canonical_genotypes.parquet").write_bytes(b"tampered synthetic bytes")
    descriptor = tmp_path / "fabricated-pedigree.json"
    descriptor.write_text(
        json.dumps(
            {
                "schema": "genome-evidence-pedigree/v1",
                "family_id": "fabricated-family",
                "members": [
                    {
                        "member_id": "fabricated-parent",
                        "subject_id": "subject-0001",
                        "m2_run": str(run),
                    },
                    {
                        "member_id": "fabricated-child",
                        "subject_id": "subject-0002",
                        "m2_run": str(child_run),
                    },
                ],
                "relationships": [
                    {
                        "assertion_id": "fabricated-edge",
                        "parent_member_id": "fabricated-parent",
                        "child_member_id": "fabricated-child",
                        "source": "user_declared",
                    }
                ],
            }
        )
    )
    try:
        analyze_family(descriptor, tmp_path / "output")
    except ValueError as error:
        assert "checksum" in str(error)
    else:
        raise AssertionError("tampered M2 input was accepted")
