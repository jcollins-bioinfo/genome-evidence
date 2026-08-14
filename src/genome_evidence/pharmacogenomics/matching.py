"""Conservative exact-identity candidate-diplotype enumeration.

The generic matcher operates only on explicitly supported autosomal diploid
small-variant definitions. It does not phase, impute, rank, or fill missing loci.
Every compatible unordered pair is retained in canonical order.
"""

import json
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any

import polars as pl

from .models import (
    CandidateDiplotype,
    ConstraintState,
    GeneInference,
    GeneOutcome,
    GeneStrategy,
    LocusEvidence,
    LocusEvidenceState,
    MatchingLimits,
)


def _id(prefix: str, value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return f"{prefix}-{sha256(payload.encode()).hexdigest()}"


def enumerate_diplotypes(
    bundle: Path,
    gene_id: str,
    evidence: tuple[LocusEvidence, ...],
    limits: MatchingLimits,
) -> GeneInference:
    """Enumerate every unordered candidate pair compatible with observed evidence.

    Parameters
    ----------
    bundle:
        Already validated bundle directory.
    gene_id:
        Explicitly selected bundle gene.
    evidence:
        Exact-locus evidence derived from M2 called observations.
    limits:
        Bounds checked before enumeration; candidates are never truncated.

    Returns
    -------
    GeneInference
        Complete stable candidate set or a fail-closed outcome.

    Notes
    -----
    Runtime is O(A²L), where A is active allele count and L is modeled locus
    count. Missing evidence can retain candidates but can never resolve one or
    create a reference allele. No score or population frequency is used.
    """
    genes = {row["gene_id"]: row for row in pl.read_parquet(bundle / "genes.parquet").to_dicts()}
    if (
        gene_id not in genes
        or GeneStrategy(genes[gene_id]["strategy"]) != GeneStrategy.STAR_HAPLOTYPE_SMALL_VARIANT
    ):
        return GeneInference(
            gene_id=gene_id, outcome=GeneOutcome.UNSUPPORTED_GENE_OR_METHOD, locus_evidence=evidence
        )
    alleles = sorted(
        row["allele_id"]
        for row in pl.read_parquet(bundle / "alleles.parquet").to_dicts()
        if row["gene_id"] == gene_id and row["status"] == "active"
    )
    loci = {
        row["locus_id"]: row
        for row in pl.read_parquet(bundle / "loci.parquet").to_dicts()
        if row["gene_id"] == gene_id
    }
    pair_count = len(alleles) * (len(alleles) + 1) // 2
    if (
        len(alleles) > limits.max_alleles
        or len(loci) > limits.max_loci
        or pair_count > limits.max_candidate_pairs
    ):
        return GeneInference(
            gene_id=gene_id,
            outcome=GeneOutcome.COMBINATORIAL_LIMIT_EXCEEDED,
            locus_evidence=evidence,
        )
    constraints: dict[tuple[str, str], ConstraintState] = {}
    for row in pl.read_parquet(bundle / "allele_locus_constraints.parquet").to_dicts():
        if row["allele_id"] in alleles and row["locus_id"] in loci:
            constraints[(row["allele_id"], row["locus_id"])] = ConstraintState(row["state"])
    observed = {item.locus_id: item for item in evidence}
    conflict_states = {
        LocusEvidenceState.DUPLICATE_CONFLICT,
        LocusEvidenceState.ALLELE_INCOMPATIBLE,
        LocusEvidenceState.PLOIDY_UNSUPPORTED,
    }
    if any(item.state in conflict_states for item in evidence):
        return GeneInference(
            gene_id=gene_id, outcome=GeneOutcome.CONFLICTING_OBSERVATIONS, locus_evidence=evidence
        )
    if any(item.state == LocusEvidenceState.UNMODELED_OBSERVED_VARIANT for item in evidence):
        return GeneInference(
            gene_id=gene_id,
            outcome=GeneOutcome.UNMODELED_OBSERVED_VARIATION,
            locus_evidence=evidence,
        )
    candidates: list[CandidateDiplotype] = []
    for first, second in combinations_with_replacement(alleles, 2):
        contradicted = False
        missing = False
        for locus_id, locus in loci.items():
            states = (constraints.get((first, locus_id)), constraints.get((second, locus_id)))
            if None in states or ConstraintState.NOT_CONSTRAINING in states:
                continue
            expected_alt = sum(state == ConstraintState.REQUIRED_ALTERNATE for state in states)
            item = observed.get(locus_id)
            if item is None or item.state in {
                LocusEvidenceState.MISSING_OR_UNASSAYED,
                LocusEvidenceState.NO_CALL,
                LocusEvidenceState.MAPPING_UNRESOLVED,
            }:
                missing = True
                continue
            actual_alt = sum(allele == locus["alternate"] for allele in item.alleles)
            if len(item.alleles) != 2 or actual_alt != expected_alt:
                contradicted = True
                break
        if not contradicted:
            candidates.append(
                CandidateDiplotype(
                    candidate_id=_id("pgx-candidate", [gene_id, first, second]),
                    gene_id=gene_id,
                    allele_a=first,
                    allele_b=second,
                    fully_evaluated=not missing,
                )
            )
    if not candidates:
        outcome = GeneOutcome.NO_COMPATIBLE_DEFINITION
    elif any(not candidate.fully_evaluated for candidate in candidates):
        outcome = GeneOutcome.INSUFFICIENT_COVERAGE
    elif len(candidates) == 1:
        outcome = GeneOutcome.RESOLVED_CANDIDATE
    else:
        outcome = GeneOutcome.AMBIGUOUS_CANDIDATES
    return GeneInference(
        gene_id=gene_id, outcome=outcome, candidates=tuple(candidates), locus_evidence=evidence
    )


def observed_locus_evidence(
    bundle: Path, variants: list[dict[str, Any]], genotypes: list[dict[str, Any]], gene_id: str
) -> tuple[LocusEvidence, ...]:
    """Align M2 called genotypes by exact assembly/chromosome/position/REF/ALT.

    Missing rows remain ``missing_or_unassayed``. Duplicate observations retain
    every genotype and observation lineage reference and fail closed on conflict.
    """
    loci = [
        row
        for row in pl.read_parquet(bundle / "loci.parquet").to_dicts()
        if row["gene_id"] == gene_id
    ]
    variant_by_id = {row["variant_id"]: row for row in variants}
    by_key: dict[tuple[str, str, int, str, str], list[dict[str, Any]]] = {}
    for genotype in genotypes:
        variant = variant_by_id.get(genotype["variant_id"])
        if variant:
            key = (
                variant["assembly"],
                variant["chromosome"],
                variant["position"],
                variant["reference"],
                variant["alternate"],
            )
            by_key.setdefault(key, []).append(genotype)
    result: list[LocusEvidence] = []
    for locus in loci:
        key = (
            locus["assembly"],
            locus["chromosome"],
            locus["position"],
            locus["reference"],
            locus["alternate"],
        )
        rows = by_key.get(key, [])
        if not rows:
            result.append(
                LocusEvidence(
                    locus_id=locus["locus_id"], state=LocusEvidenceState.MISSING_OR_UNASSAYED
                )
            )
            continue
        allele_sets = {tuple(sorted(row["alleles"])) for row in rows}
        alleles = next(iter(allele_sets))
        if len(allele_sets) > 1:
            state = LocusEvidenceState.DUPLICATE_CONFLICT
        elif len(rows) > 1:
            state = LocusEvidenceState.DUPLICATE_CONCORDANT
        elif len(alleles) != 2:
            state = LocusEvidenceState.PLOIDY_UNSUPPORTED
        else:
            alt_count = sum(value == locus["alternate"] for value in alleles)
            state = (
                LocusEvidenceState.OBSERVED_REFERENCE,
                LocusEvidenceState.OBSERVED_HETEROZYGOUS,
                LocusEvidenceState.OBSERVED_ALTERNATE,
            )[alt_count]
        result.append(
            LocusEvidence(
                locus_id=locus["locus_id"],
                state=state,
                alleles=alleles,
                genotype_ids=tuple(sorted(row["genotype_id"] for row in rows)),
                observation_references=tuple(sorted(row["observation_reference"] for row in rows)),
            )
        )
    return tuple(result)
