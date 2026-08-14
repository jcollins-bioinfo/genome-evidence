import ast
import json
import subprocess
import sys
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
