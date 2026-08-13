"""Explicit local, checksummed M2 reference providers."""

import json
from hashlib import sha256
from pathlib import Path
from typing import Protocol

from pydantic import BaseModel, ConfigDict, Field

from genome_evidence.normalization.models import ResourceIdentity, StrandTransform


def canonical_assembly(value: str) -> str | None:
    return {"grch37": "GRCh37", "hg19": "GRCh37", "grch38": "GRCh38", "hg38": "GRCh38"}.get(
        value.lower()
    )


def canonical_chromosome(value: str) -> str:
    token = value.removeprefix("chr").upper()
    return "MT" if token == "M" else token


class MarkerDefinition(BaseModel):
    model_config = ConfigDict(frozen=True)
    marker_id: str
    assembly: str
    chromosome: str
    position: int = Field(gt=0)
    reference: str
    alternate: str
    orientation: StrandTransform
    orientation_authoritative: bool = False
    rsid: str | None = None


class MarkerProvider(Protocol):
    identity: ResourceIdentity

    def definitions(self, marker_id: str) -> tuple[MarkerDefinition, ...]: ...


class ReferenceProvider(Protocol):
    identity: ResourceIdentity

    def sequence(self, chromosome: str, position: int, length: int) -> str | None: ...


class LiftoverProvider(Protocol):
    identity: ResourceIdentity

    def lift(self, chromosome: str, zero_based_position: int) -> tuple[tuple[str, int], ...]: ...


def _identity(path: Path, kind: str, name: str, version: str, **kwargs: str) -> ResourceIdentity:
    return ResourceIdentity(
        resource_type=kind,
        logical_name=name,
        version=version,
        sha256=sha256(path.read_bytes()).hexdigest(),
        local_identity=path.name,
        **kwargs,
    )


class JsonMarkerProvider:
    def __init__(self, path: Path, name: str = "marker-definitions", version: str = "1") -> None:
        self.identity = _identity(path, "marker_definitions", name, version)
        rows = json.loads(path.read_text())
        self._rows = tuple(MarkerDefinition.model_validate(row) for row in rows)

    def definitions(self, marker_id: str) -> tuple[MarkerDefinition, ...]:
        return tuple(row for row in self._rows if row.marker_id == marker_id)


class FastaReferenceProvider:
    def __init__(self, path: Path, assembly: str, version: str = "1") -> None:
        self.identity = _identity(path, "reference_sequence", path.stem, version, assembly=assembly)
        sequences: dict[str, str] = {}
        current: str | None = None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                current = canonical_chromosome(line[1:].split()[0])
                sequences[current] = ""
            elif current:
                sequences[current] += line.strip().upper()
        self._sequences = sequences

    def sequence(self, chromosome: str, position: int, length: int) -> str | None:
        seq = self._sequences.get(canonical_chromosome(chromosome))
        if seq is None or position < 1 or position - 1 + length > len(seq):
            return None
        return seq[position - 1 : position - 1 + length]


class JsonLiftoverProvider:
    """Tiny local mapping format; input/output positions are explicitly zero-based."""

    def __init__(self, path: Path, source: str, target: str, version: str = "1") -> None:
        self.identity = _identity(
            path, "liftover", path.stem, version, source_assembly=source, target_assembly=target
        )
        self._rows = json.loads(path.read_text())

    def lift(self, chromosome: str, zero_based_position: int) -> tuple[tuple[str, int], ...]:
        key = f"{canonical_chromosome(chromosome)}:{zero_based_position}"
        return tuple((canonical_chromosome(x[0]), int(x[1])) for x in self._rows.get(key, []))
