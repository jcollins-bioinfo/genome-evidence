import json
from hashlib import sha256
from pathlib import Path

from genome_evidence.normalization.resources import FastaReferenceProvider, JsonMarkerProvider


def test_indexed_fasta_reads_across_wrapped_lines_and_records_index(tmp_path: Path) -> None:
    fasta = tmp_path / "GRCh38.fa"
    fasta.write_bytes(b">chr1\nACGT\nTGCA\n")
    index = tmp_path / "GRCh38.fa.fai"
    index.write_text("chr1\t8\t6\t4\t5\n")

    provider = FastaReferenceProvider(fasta, "GRCh38", "synthetic-indexed")

    assert provider.sequence("1", 3, 4) == "GTTG"
    assert provider.sequence("1", 8, 2) is None
    assert provider.identity.index_sha256 == sha256(index.read_bytes()).hexdigest()
    provider.close()


def test_marker_provider_indexes_duplicate_definitions(tmp_path: Path) -> None:
    path = tmp_path / "markers.json"
    row = {
        "marker_id": "synthetic",
        "assembly": "GRCh38",
        "chromosome": "1",
        "position": 1,
        "reference": "A",
        "alternate": "C",
        "orientation": "none",
        "orientation_authoritative": True,
    }
    path.write_text(json.dumps([row, {**row, "position": 2}]))

    provider = JsonMarkerProvider(path)

    assert [definition.position for definition in provider.definitions("synthetic")] == [1, 2]
    assert provider.definitions("absent") == ()
