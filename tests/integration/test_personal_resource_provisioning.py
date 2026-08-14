import gzip
import json
import subprocess
from hashlib import md5
from pathlib import Path

import pytest

import genome_evidence.workspace.resources as resource_module
from genome_evidence.workspace import (
    WorkspaceConfig,
    import_23andme_source,
    initialize_workspace,
    provision_personal_normalization_resources,
    run_personal_m1_m2,
)


def test_provisioned_resources_drive_personal_m1_m2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = initialize_workspace(tmp_path / "private", WorkspaceConfig(subject_id="subject-0001"))
    inbox = root / "inputs/inbox/23andme/source.txt"
    inbox.write_text("# genome build: GRCh37\nrs1\t1\t2\tCG\nrsMissing\t1\t8\tAA\n")
    import_23andme_source(root, inbox, "subject-0001")
    fasta_payload = gzip.compress(b">chr1\nAAACAAAAAA\n")
    monkeypatch.setattr(
        resource_module,
        "FASTA_UPSTREAM_MD5",
        md5(fasta_payload, usedforsecurity=False).hexdigest(),
    )

    def fake_download(url: str, destination: Path) -> None:
        if url == resource_module.KENT_TOOL_URL:
            destination.write_bytes(b"synthetic kent executable")
        elif url == resource_module.FASTA_URL:
            destination.write_bytes(fasta_payload)
        else:
            raise AssertionError(f"unexpected full download: {url}")

    def fake_runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if len(args) == 1:
            return subprocess.CompletedProcess(
                args, 255, stdout="", stderr="### kent source version 479 ###\n"
            )
        output = Path(args[-1])
        if args[2] == resource_module.DBSNP_URLS["GRCh37"]:
            output.write_text("chr1\t1\t2\trs1\tC\t1\tG,\t0\n")
        elif args[2] == resource_module.DBSNP_URLS["GRCh38"]:
            output.write_text("chr1\t3\t4\trs1\tC\t1\tG,\t0\n")
        else:
            raise AssertionError(f"unexpected indexed query: {args[2]}")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    provisioned = provision_personal_normalization_resources(
        root,
        {},
        working_root=tmp_path / "work",
        download=fake_download,
        runner=fake_runner,
    )

    assert provisioned.marker_count == 1
    assert provisioned.mapped_marker_count == 1
    assert provisioned.unresolved_marker_count == 1
    selection = json.loads((root / "config/normalization_resources.json").read_text())
    assert selection["source_assembly"] == "GRCh37"
    assert selection["marker_definitions"].startswith("references/markers/23andme/")
    assert selection["grch38_fasta"].startswith("references/genome/grch38/")
    assert selection["grch37_to_grch38_liftover"].startswith(
        "references/liftover/grch37_to_grch38/"
    )

    normalized = run_personal_m1_m2(
        root, "subject-0001", {}, working_root=tmp_path / "normalization-work"
    )

    assert normalized.observation_count == 2
    assert normalized.mapping_count == 2
    assert normalized.canonical_genotype_count == 1
