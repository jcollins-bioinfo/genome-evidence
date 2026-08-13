"""Source-agnostic, local-only phasing reference validation."""

from hashlib import sha256
from pathlib import Path

from pydantic import ValidationError

from .models import ReferenceManifest, ReferenceValidation


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_phasing_reference(bundle_directory: Path) -> ReferenceValidation:
    root = bundle_directory.resolve()
    path = root / "manifest.json"
    try:
        manifest = ReferenceManifest.model_validate_json(path.read_bytes())
    except (OSError, ValidationError) as error:
        raise ValueError("invalid phasing reference manifest") from error
    for chromosome in manifest.per_chromosome:
        declared = {
            chromosome.panel,
            chromosome.panel_index,
            chromosome.genetic_map,
            chromosome.variants,
        }
        if set(chromosome.artifacts) != declared:
            raise ValueError(f"chromosome {chromosome.chromosome} artifact inventory mismatch")
        for relative, identity in chromosome.artifacts.items():
            artifact = (root / relative).resolve()
            if not artifact.is_relative_to(root) or not artifact.is_file():
                raise ValueError("reference artifact path is missing or escapes bundle")
            if artifact.stat().st_size != identity.byte_size or _hash(artifact) != identity.sha256:
                raise ValueError(f"reference artifact integrity failure: {relative}")
        rows = [
            line.split("\t")
            for line in (root / chromosome.variants).read_text().splitlines()
            if line and not line.startswith("#")
        ]
        keys = []
        for row in rows:
            if len(row) != 6 or row[0] != chromosome.chromosome:
                raise ValueError("invalid canonical variant index")
            pos, ref, alt, ac, an = int(row[1]), row[2], row[3], int(row[4]), int(row[5])
            if pos < 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt or not 0 <= ac <= an:
                raise ValueError("invalid canonical reference allele/count")
            keys.append((pos, ref, alt))
        if keys != sorted(set(keys)):
            raise ValueError("canonical reference variants must be sorted and unique")
        maps = [
            line.split()
            for line in (root / chromosome.genetic_map).read_text().splitlines()
            if line and not line.startswith("#")
        ]
        points = [(int(row[0]), float(row[1])) for row in maps]
        if (
            len(points) < 2
            or points != sorted(points)
            or any(b[1] < a[1] for a, b in zip(points, points[1:], strict=False))
        ):
            raise ValueError("genetic map must have monotonic position and cM coverage")
    return ReferenceValidation(directory=root, manifest=manifest, manifest_sha256=_hash(path))
