import gzip
import json
import subprocess
from hashlib import md5
from pathlib import Path

import pytest

import genome_evidence.workspace.resources as resource_module
from genome_evidence.workspace import (
    ProvisioningIncomplete,
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
            raise OSError(f"synthetic unavailable resource: {url}")

    def fake_runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if len(args) == 1:
            return subprocess.CompletedProcess(
                args,
                255,
                stdout="",
                stderr=(
                    "bigBedNamedItems - Extract item of given name from bigBed\n"
                    "usage:\n   bigBedNamedItems file.bb name output.bed\n"
                    "   -nameFile - treat name as a file\n"
                ),
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
        {"GENOME_EVIDENCE_DBSNP_COVERAGE_POLICY": "full_remote_v1"},
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

    fasta_index = root / f"{provisioned.selection.grch38_fasta}.fai"
    fasta_index.unlink()
    with pytest.raises(ValueError, match="failed validation|missing or unsafe"):
        resource_module.load_normalization_resource_selection(
            root,
            provisioned.selection.source_sha256,
        )


def test_provisioning_resumes_after_target_dbsnp_interruption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = initialize_workspace(tmp_path / "private", WorkspaceConfig(subject_id="subject-0001"))
    private_rsid = "rs987654321"
    private_genotype = "CG"
    inbox = root / "inputs/inbox/23andme/source.txt"
    inbox.write_text(f"# genome build: GRCh37\n{private_rsid}\t1\t2\t{private_genotype}\n")
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
            raise OSError(f"synthetic unavailable resource: {url}")

    query_counts = {"source": 0, "target": 0}

    def probe(args: list[str]) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args,
            255,
            stdout="",
            stderr=(
                "bigBedNamedItems - Extract item of given name from bigBed\n"
                "usage:\n   bigBedNamedItems file.bb name output.bed\n"
                "   -nameFile - treat name as a file\n"
            ),
        )

    def interrupted_runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if len(args) == 1:
            return probe(args)
        output = Path(args[-1])
        if args[2] == resource_module.DBSNP_URLS["GRCh37"]:
            query_counts["source"] += 1
            output.write_text(f"chr1\t1\t2\t{private_rsid}\tC\t1\tG,\t0\n")
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")
        if args[2] == resource_module.DBSNP_URLS["GRCh38"]:
            query_counts["target"] += 1
            return subprocess.CompletedProcess(
                args,
                255,
                stdout="",
                stderr="transient remote indexed-query failure",
            )
        raise AssertionError(f"unexpected indexed query: {args[2]}")

    environment = {
        "GENOME_EVIDENCE_DBSNP_COVERAGE_POLICY": "full_remote_v1",
        "GENOME_EVIDENCE_DBSNP_BATCH_SIZE": "250",
        "GENOME_EVIDENCE_QUERY_ATTEMPTS": "1",
        "GENOME_EVIDENCE_QUERY_TIMEOUT_SECONDS": "60",
    }
    with pytest.raises(ProvisioningIncomplete) as failure:
        provision_personal_normalization_resources(
            root,
            environment,
            working_root=tmp_path / "work",
            download=fake_download,
            runner=interrupted_runner,
            sleep=lambda _seconds: None,
        )

    assert query_counts == {"source": 1, "target": 1}
    assert not (root / "config/normalization_resources.json").exists()
    assert failure.value.checkpoint_path.is_dir()
    assert failure.value.log_path.is_file()

    def resumed_runner(args: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        if len(args) == 1:
            return probe(args)
        output = Path(args[-1])
        if args[2] == resource_module.DBSNP_URLS["GRCh37"]:
            raise AssertionError("completed source query checkpoint must be reused")
        if args[2] == resource_module.DBSNP_URLS["GRCh38"]:
            query_counts["target"] += 1
            output.write_text(f"chr1\t3\t4\t{private_rsid}\tC\t1\tG,\t0\n")
            return subprocess.CompletedProcess(args, 0, stdout="", stderr="")
        raise AssertionError(f"unexpected indexed query: {args[2]}")

    provisioned = provision_personal_normalization_resources(
        root,
        environment,
        working_root=tmp_path / "work",
        download=fake_download,
        runner=resumed_runner,
        sleep=lambda _seconds: None,
    )

    assert query_counts == {"source": 1, "target": 2}
    assert provisioned.selection_path.is_file()
    assert provisioned.selection.source_assembly == "GRCh37"
    log_text = provisioned.log_path.read_text(encoding="utf-8")
    assert '"event": "dbsnp.grch37.resume"' in log_text
    assert '"event": "session.complete"' in log_text
    assert private_rsid not in log_text
    assert private_genotype not in log_text

    selector = root / "config/normalization_resources.json"
    selector_completion = root / "config/normalization_resources.COMPLETED.json"
    selector_publishing = root / "config/normalization_resources.PUBLISHING.json"
    selector_completion.unlink()
    selector.write_text("{", encoding="utf-8")
    selector_publishing.write_text("{}\n", encoding="utf-8")

    def checkpoint_only_runner(
        args: list[str], **_kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        if len(args) == 1:
            return probe(args)
        raise AssertionError("all completed dbSNP query checkpoints must be reused")

    recovered = provision_personal_normalization_resources(
        root,
        environment,
        working_root=tmp_path / "work",
        download=fake_download,
        runner=checkpoint_only_runner,
        sleep=lambda _seconds: None,
    )

    assert recovered.selection_path == selector
    assert selector_completion.is_file()
    assert not selector_publishing.exists()
    assert json.loads(selector.read_text())["source_sha256"] == recovered.selection.source_sha256

    provenance_relative = Path(recovered.selection.provenance_manifest)
    bundle_completion = root / provenance_relative.with_name(
        f"{provenance_relative.stem}.COMPLETED.json"
    )
    marker_definitions = root / recovered.selection.marker_definitions
    selector.unlink()
    selector_completion.unlink()
    bundle_completion.unlink()
    marker_definitions.write_text("{", encoding="utf-8")

    repaired = provision_personal_normalization_resources(
        root,
        environment,
        working_root=tmp_path / "work",
        download=fake_download,
        runner=checkpoint_only_runner,
        sleep=lambda _seconds: None,
    )

    assert repaired.selection_path.is_file()
    assert bundle_completion.is_file()
    assert isinstance(json.loads(marker_definitions.read_text()), list)
