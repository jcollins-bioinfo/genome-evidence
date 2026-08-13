import json
from pathlib import Path


def test_notebook_portfolio_and_profiles_are_synchronized() -> None:
    expected = [
        "00_initialize_private_workspace.ipynb",
        "01_ingest_and_normalize_genome.ipynb",
        "02_versioned_external_evidence.ipynb",
        "03_evidence_oriented_clinical_review.ipynb",
        "04_population_structure_projection.ipynb",
        "05_reference_panel_phasing_and_imputation.ipynb",
    ]
    assert [p.name for p in sorted(Path("notebooks").glob("*.ipynb"))] == expected
    root, index = Path("README.md").read_text(), Path("notebooks/README.md").read_text()
    for name in expected:
        assert root.count(f"notebooks/{name}") == 2
        assert index.count(name) == 3
        notebook = json.loads((Path("notebooks") / name).read_text())
        source = "".join("".join(c.get("source", [])) for c in notebook["cells"])
        assert "GENOME_EVIDENCE_PROFILE" in source
        assert "personal_drive" in source and "synthetic_ci" in source
        assert all(
            c.get("execution_count") is None for c in notebook["cells"] if c["cell_type"] == "code"
        )
        assert all(not c.get("outputs") for c in notebook["cells"] if c["cell_type"] == "code")
    assert "01_ingest_and_normalize_synthetic_genome" not in root + index
