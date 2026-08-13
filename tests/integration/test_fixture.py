from pathlib import Path

import polars as pl


def test_synthetic_fixture_loads() -> None:
    fixture = Path(__file__).parents[1] / "fixtures" / "synthetic_genotypes.tsv"
    frame = pl.read_csv(fixture, separator="\t", comment_prefix="#")
    assert frame.height == 3
    missing = frame.filter(pl.col("marker_id") == "synthetic_marker_missing")
    assert missing.item(0, "genotype") == "--"
