import ast
import json
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Any, NoReturn

import pytest


def canonical_bootstrap() -> str:
    tree = ast.parse(Path("scripts/sync_notebook_bootstrap.py").read_text())
    assignment = next(
        node
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "BOOTSTRAP" for target in node.targets
        )
    )
    assert isinstance(assignment, ast.Assign)
    value = ast.literal_eval(assignment.value)
    assert isinstance(value, str)
    return value


def test_notebook_portfolio_and_profiles_are_synchronized() -> None:
    expected = [
        "00_initialize_private_workspace.ipynb",
        "01_ingest_and_normalize_genome.ipynb",
        "02_versioned_external_evidence.ipynb",
        "03_evidence_oriented_clinical_review.ipynb",
        "04_population_structure_projection.ipynb",
        "05_reference_panel_phasing_and_imputation.ipynb",
        "06_polygenic_score_calculation.ipynb",
        "07_pharmacogenomics_evidence.ipynb",
    ]
    assert [p.name for p in sorted(Path("notebooks").glob("*.ipynb"))] == expected
    root, index = Path("README.md").read_text(), Path("notebooks/README.md").read_text()
    for name in expected:
        assert root.count(f"notebooks/{name}") == 2
        assert index.count(name) == 3
        notebook = json.loads((Path("notebooks") / name).read_text())
        tags = [c.get("metadata", {}).get("tags", []) for c in notebook["cells"]]
        assert [
            next(i for i, x in enumerate(tags) if tag in x)
            for tag in ("parameters", "genome-evidence-bootstrap", "genome-evidence-workspace")
        ] == [1, 2, 3]
        source = "".join("".join(c.get("source", [])) for c in notebook["cells"])
        assert "GENOME_EVIDENCE_PROFILE" in source
        assert "personal_drive" in source and "synthetic_ci" in source
        assert all(
            c.get("execution_count") is None for c in notebook["cells"] if c["cell_type"] == "code"
        )
        assert all(not c.get("outputs") for c in notebook["cells"] if c["cell_type"] == "code")
    assert "01_ingest_and_normalize_synthetic_genome" not in root + index


def test_bootstrap_cells_are_byte_identical() -> None:
    notebooks = [json.loads(p.read_text()) for p in sorted(Path("notebooks").glob("*.ipynb"))]
    cells = [
        "".join(
            next(
                c
                for c in n["cells"]
                if "genome-evidence-bootstrap" in c.get("metadata", {}).get("tags", [])
            )["source"]
        )
        for n in notebooks
    ]
    assert len(set(cells)) == 1
    assert cells[0] == canonical_bootstrap()


def test_synthetic_bootstrap_uses_installed_package_without_external_operations(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def forbidden(*args: object, **kwargs: object) -> NoReturn:
        raise AssertionError("synthetic bootstrap attempted an external process")

    monkeypatch.setattr(subprocess, "run", forbidden)
    namespace: dict[str, Any] = {
        "PROFILE": "synthetic_ci",
        "REPOSITORY_URL": "https://github.com/jcollins-bioinfo/genome-evidence.git",
        "REPOSITORY_REF": "main",
        "WORKSPACE_ROOT": "/must/not/be/touched",
        "SUBJECT_ID": "subject-0001",
        "sys": sys,
    }
    exec(canonical_bootstrap(), namespace)
    assert namespace["BOOTSTRAP_STATUS"]["profile"] == "synthetic_ci"
    assert namespace["BOOTSTRAP_STATUS"]["import_path"] == "installed-ci-package"


def test_personal_bootstrap_resolves_fetched_head_not_a_local_branch() -> None:
    source = canonical_bootstrap()
    assert '["git", "-C", str(CHECKOUT), "fetch", "--force", "origin", REPOSITORY_REF]' in source
    assert '"FETCH_HEAD^{commit}"' in source
    assert 'f"{REPOSITORY_REF}^{{commit}}"' not in source
    assert source.index('importlib.import_module("genome_evidence")') > source.index(
        '"pip",\n            "install"'
    )


def test_personal_notebook_workflow_publishes_and_resolves_runs() -> None:
    notebook_01 = Path("notebooks/01_ingest_and_normalize_genome.ipynb").read_text()
    notebook_04 = Path("notebooks/04_population_structure_projection.ipynb").read_text()

    assert "run_personal_m1_m2" in notebook_01
    assert "GENOME_EVIDENCE_NORMALIZATION_RUN" in notebook_01
    assert "resolve_personal_m2_run" in notebook_04
    assert "resolve_personal_population_bundle" in notebook_04
    assert "publish_completed_run" in notebook_04
    assert "M2 run missing; complete notebook 01" not in notebook_04


def test_personal_bootstrap_imports_checkout_in_same_fresh_process(tmp_path: Path) -> None:
    checkout = tmp_path / "checkout"
    package = checkout / "src/genome_evidence"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text('__version__ = "0.1.0"\n')
    metadata = checkout / "src/genome_evidence-0.1.0.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: genome-evidence\nVersion: 0.1.0\n"
    )
    (checkout / "uv.lock").write_text("synthetic bootstrap test lock\n")

    source = canonical_bootstrap().replace(
        'CHECKOUT = Path("/content/genome-evidence-src")',
        f"CHECKOUT = Path({str(checkout)!r})",
    )
    harness = tmp_path / "fresh_personal_bootstrap.py"
    harness.write_text(
        textwrap.dedent(
            f"""
            import json
            import subprocess
            import sys

            repository_url = "https://github.com/jcollins-bioinfo/genome-evidence.git"
            resolved_commit = "0123456789abcdef0123456789abcdef01234567"
            calls = []

            def fake_run(args, **kwargs):
                calls.append(args)
                if args[-3:] == ["remote", "get-url", "origin"]:
                    stdout = repository_url + "\\n"
                elif args[-2:] == ["status", "--porcelain"]:
                    stdout = ""
                elif args[-3:] == ["rev-parse", "--verify", "FETCH_HEAD^{{commit}}"]:
                    stdout = resolved_commit + "\\n"
                else:
                    stdout = ""
                return subprocess.CompletedProcess(args, 0, stdout=stdout)

            subprocess.run = fake_run
            namespace = {{
                "PROFILE": "personal_drive",
                "REPOSITORY_URL": repository_url,
                "REPOSITORY_REF": "main",
                "WORKSPACE_ROOT": "/not-used",
                "SUBJECT_ID": "subject-0001",
                "sys": sys,
            }}
            exec({source!r}, namespace)
            status = namespace["BOOTSTRAP_STATUS"]
            assert status["resolved_commit"] == resolved_commit
            assert status["import_path"] == "src/genome_evidence/__init__.py"
            assert any("pip" in call and "install" in call for call in calls)
            print(json.dumps(status, sort_keys=True))
            """
        )
    )

    completed = subprocess.run(
        [sys.executable, str(harness)],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    status = json.loads(completed.stdout.splitlines()[-1])
    assert status["profile"] == "personal_drive"
