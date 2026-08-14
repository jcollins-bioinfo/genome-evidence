"""Explicit local, checksummed M2 reference providers."""

import json
import os
from hashlib import sha256
from pathlib import Path
from typing import BinaryIO, Protocol

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


def _hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _identity(
    path: Path,
    kind: str,
    name: str,
    version: str,
    *,
    index_sha256: str | None = None,
    **kwargs: str,
) -> ResourceIdentity:
    return ResourceIdentity(
        resource_type=kind,
        logical_name=name,
        version=version,
        sha256=_hash(path),
        index_sha256=index_sha256,
        local_identity=path.name,
        **kwargs,
    )


class JsonMarkerProvider:
    def __init__(self, path: Path, name: str = "marker-definitions", version: str = "1") -> None:
        self.identity = _identity(path, "marker_definitions", name, version)
        rows = json.loads(path.read_text())
        indexed: dict[str, list[MarkerDefinition]] = {}
        for row in rows:
            definition = MarkerDefinition.model_validate(row)
            indexed.setdefault(definition.marker_id, []).append(definition)
        self._rows = {key: tuple(value) for key, value in indexed.items()}

    def definitions(self, marker_id: str) -> tuple[MarkerDefinition, ...]:
        return self._rows.get(marker_id, ())


class FastaReferenceProvider:
    def __init__(self, path: Path, assembly: str, version: str = "1") -> None:
        index_path = Path(f"{path}.fai")
        index_hash = _hash(index_path) if index_path.is_file() else None
        self.identity = _identity(
            path,
            "reference_sequence",
            path.stem,
            version,
            index_sha256=index_hash,
            assembly=assembly,
        )
        self._handle: BinaryIO | None = None
        self._index: dict[str, tuple[int, int, int, int]] = {}
        sequences: dict[str, str] = {}
        if index_hash is not None:
            for line in index_path.read_text().splitlines():
                columns = line.split("\t")
                if len(columns) < 5:
                    raise ValueError("invalid FASTA index row")
                chromosome = canonical_chromosome(columns[0])
                if chromosome in self._index:
                    raise ValueError("duplicate chromosome in FASTA index")
                try:
                    length, offset, line_bases, line_width = map(int, columns[1:5])
                except ValueError as error:
                    raise ValueError("invalid FASTA index coordinates") from error
                if min(length, offset) < 0 or line_bases <= 0 or line_width < line_bases:
                    raise ValueError("invalid FASTA index dimensions")
                self._index[chromosome] = (length, offset, line_bases, line_width)
            self._handle = path.open("rb")
        else:
            if path.stat().st_size > 50 * 1024 * 1024:
                raise ValueError("large FASTA requires a matching .fai index")
            current: str | None = None
            fragments: dict[str, list[str]] = {}
            for line in path.read_text().splitlines():
                if line.startswith(">"):
                    current = canonical_chromosome(line[1:].split()[0])
                    fragments[current] = []
                elif current:
                    fragments[current].append(line.strip().upper())
            sequences = {name: "".join(parts) for name, parts in fragments.items()}
        self._sequences = sequences

    def sequence(self, chromosome: str, position: int, length: int) -> str | None:
        chromosome = canonical_chromosome(chromosome)
        if self._handle is not None:
            entry = self._index.get(chromosome)
            if entry is None or position < 1 or length < 0 or position - 1 + length > entry[0]:
                return None
            _, offset, line_bases, line_width = entry
            cursor = position - 1
            remaining = length
            chunks = []
            while remaining:
                within_line = cursor % line_bases
                count = min(remaining, line_bases - within_line)
                file_offset = offset + (cursor // line_bases) * line_width + within_line
                chunk = os.pread(self._handle.fileno(), count, file_offset)
                if len(chunk) != count:
                    return None
                chunks.append(chunk)
                cursor += count
                remaining -= count
            try:
                return b"".join(chunks).decode("ascii").upper()
            except UnicodeDecodeError:
                return None
        seq = self._sequences.get(canonical_chromosome(chromosome))
        if seq is None or position < 1 or position - 1 + length > len(seq):
            return None
        return seq[position - 1 : position - 1 + length]

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def __del__(self) -> None:
        handle = getattr(self, "_handle", None)
        if handle is not None:
            handle.close()


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
