"""Pure autosomal diploid Mendelian transmission enumeration for M9."""

from hashlib import sha256

from .models import (
    CompatibilityStatus,
    Informativeness,
    SegregationEvidence,
    TransmissionAssignment,
    TransmissionStatus,
)


def evaluate_site(
    *,
    variant_id: str,
    child_member_id: str,
    child_genotype: tuple[str, str] | None,
    parent_genotypes: tuple[tuple[str, str] | None, ...],
    relationship_assertion_ids: tuple[str, ...],
    genotype_record_ids: tuple[str, ...],
    reference: str,
    alternate: str,
    chromosome: str,
    conflicting: bool = False,
    unsupported_ploidy: bool = False,
) -> SegregationEvidence:
    """Enumerate compatible transmissions for one trio or duo site.

    The unknown contributor in a duo ranges over both supported alleles. Results
    are conditional compatibility statements only; inconsistency is not labelled
    as mutation, sample identity error, or non-parentage.
    """
    identity = (
        "m9e-"
        + sha256(
            repr(
                (
                    variant_id,
                    child_member_id,
                    relationship_assertion_ids,
                    genotype_record_ids,
                    child_genotype,
                    parent_genotypes,
                    reference,
                    alternate,
                    chromosome,
                    conflicting,
                    unsupported_ploidy,
                )
            ).encode()
        ).hexdigest()
    )
    base = dict(
        evidence_id=identity,
        variant_id=variant_id,
        child_member_id=child_member_id,
        relationship_assertion_ids=relationship_assertion_ids,
        genotype_record_ids=genotype_record_ids,
    )
    unresolved: CompatibilityStatus | None = None
    if chromosome not in {str(value) for value in range(1, 23)}:
        unresolved = CompatibilityStatus.UNSUPPORTED_LOCUS
    elif unsupported_ploidy:
        unresolved = CompatibilityStatus.UNSUPPORTED_PLOIDY
    elif len(reference) == 0 or len(alternate) == 0 or reference == alternate:
        unresolved = CompatibilityStatus.UNSUPPORTED_GENOTYPE_REPRESENTATION
    elif conflicting:
        unresolved = CompatibilityStatus.INDETERMINATE_CONFLICTING
    elif child_genotype is None or not parent_genotypes or any(x is None for x in parent_genotypes):
        unresolved = CompatibilityStatus.INDETERMINATE_MISSING
    if unresolved is not None:
        return SegregationEvidence(
            **base,
            compatibility=unresolved,
            informativeness=Informativeness.NOT_EVALUATED,
            transmission_status=(
                TransmissionStatus.UNRESOLVED_INPUT
                if unresolved
                in {
                    CompatibilityStatus.INDETERMINATE_MISSING,
                    CompatibilityStatus.INDETERMINATE_CONFLICTING,
                }
                else TransmissionStatus.NOT_APPLICABLE
            ),
        )
    alphabet = {reference, alternate}
    assert child_genotype is not None
    parents = tuple(item for item in parent_genotypes if item is not None)
    if any(set(genotype) - alphabet for genotype in (child_genotype, *parents)):
        return SegregationEvidence(
            **base,
            compatibility=CompatibilityStatus.UNSUPPORTED_GENOTYPE_REPRESENTATION,
            informativeness=Informativeness.NOT_EVALUATED,
            transmission_status=TransmissionStatus.NOT_APPLICABLE,
        )
    if len(parents) == 1:
        assignments = {
            (known, unknown)
            for known in set(parents[0])
            for unknown in alphabet
            if sorted((known, unknown)) == sorted(child_genotype)
        }
    else:
        assignments = {
            (first, second)
            for first in set(parents[0])
            for second in set(parents[1])
            if sorted((first, second)) == sorted(child_genotype)
        }
    ordered = tuple(
        TransmissionAssignment(parent_1_allele=first, parent_2_allele=second)
        for first, second in sorted(assignments)
    )
    if not ordered:
        return SegregationEvidence(
            **base,
            compatibility=CompatibilityStatus.INCONSISTENT,
            informativeness=Informativeness.INFORMATIVE,
            transmission_status=TransmissionStatus.NO_COMPATIBLE_TRANSMISSION,
        )
    unique = len(ordered) == 1
    return SegregationEvidence(
        **base,
        compatibility=CompatibilityStatus.CONSISTENT,
        informativeness=Informativeness.INFORMATIVE if unique else Informativeness.UNINFORMATIVE,
        transmission_status=(
            TransmissionStatus.UNIQUE_TRANSMISSION
            if unique
            else TransmissionStatus.AMBIGUOUS_TRANSMISSION
        ),
        assignments=ordered,
    )
