import re
from pathlib import Path


def test_every_notebook_has_exact_colab_links_in_both_indexes() -> None:
    notebooks = sorted(Path("notebooks").glob("*.ipynb"))
    root = Path("README.md").read_text()
    index = Path("notebooks/README.md").read_text()
    prefix = "https://colab.research.google.com/github/jcollins-bioinfo/genome-evidence/blob/main/"
    for notebook in notebooks:
        target = prefix + notebook.as_posix()
        assert root.count(target) == 1
        assert index.count(target) == 1
        assert notebook.is_file()
    listed = re.findall(r"\]\((\d\d_[^)]+\.ipynb)\)", index)
    assert listed == [p.name for p in notebooks]
