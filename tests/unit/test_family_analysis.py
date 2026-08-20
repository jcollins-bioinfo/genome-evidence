"""Epistemic regression tests for synthetic M9 pedigree and enumeration logic."""

from collections.abc import Callable
from itertools import product
from pathlib import Path

import pytest

from genome_evidence.family_analysis import (
    CompatibilityStatus,
    ParentRelationship,
    PedigreeDescriptor,
    PedigreeMember,
    RelationshipSource,
    TransmissionStatus,
    evaluate_site,
    validate_pedigree,
)


def _member(number: int) -> PedigreeMember:
    return PedigreeMember(
        member_id=f"fabricated-member-{number}",
        subject_id=f"subject-{number:04d}",
        m2_run=Path(f"synthetic-run-{number}"),
    )


def _edge(parent: int, child: int, suffix: str = "") -> ParentRelationship:
    return ParentRelationship(
        assertion_id=f"fabricated-edge-{parent}-{child}{suffix}",
        parent_member_id=f"fabricated-member-{parent}",
        child_member_id=f"fabricated-member-{child}",
        source=RelationshipSource.USER_DECLARED,
    )


def test_valid_trio_duo_and_order_are_deterministic() -> None:
    trio = PedigreeDescriptor(
        family_id="fabricated-family",
        members=(_member(3), _member(1), _member(2)),
        relationships=(_edge(2, 3), _edge(1, 3)),
    )
    validated = validate_pedigree(trio)
    assert [member.member_id for member in validated.members] == sorted(
        member.member_id for member in trio.members
    )
    assert (
        len(
            validate_pedigree(
                PedigreeDescriptor(
                    family_id="fabricated-duo",
                    members=(_member(1), _member(2)),
                    relationships=(_edge(1, 2),),
                )
            ).relationships
        )
        == 1
    )


@pytest.mark.parametrize(
    ("descriptor", "message"),
    [
        (
            lambda: PedigreeDescriptor(
                family_id="x", members=(_member(1), _member(1)), relationships=(_edge(1, 2),)
            ),
            "duplicate",
        ),
        (
            lambda: PedigreeDescriptor(
                family_id="x", members=(_member(1),), relationships=(_edge(1, 1),)
            ),
            "self-parent",
        ),
        (
            lambda: PedigreeDescriptor(
                family_id="x",
                members=(_member(1), _member(2)),
                relationships=(_edge(1, 2), _edge(1, 2, "-other")),
            ),
            "duplicate biological-parent",
        ),
        (
            lambda: PedigreeDescriptor(
                family_id="x", members=(_member(1), _member(2)), relationships=(_edge(1, 3),)
            ),
            "absent",
        ),
        (
            lambda: PedigreeDescriptor(
                family_id="x",
                members=(_member(1), _member(2)),
                relationships=(_edge(1, 2), _edge(2, 1)),
            ),
            "cycle",
        ),
        (
            lambda: PedigreeDescriptor(
                family_id="x",
                members=(_member(1), _member(2), _member(3), _member(4)),
                relationships=(_edge(1, 4), _edge(2, 4), _edge(3, 4)),
            ),
            "at most two",
        ),
    ],
)
def test_invalid_pedigrees_fail_closed(
    descriptor: Callable[[], PedigreeDescriptor], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        validate_pedigree(descriptor())


def test_pedigree_collapse_is_not_a_cycle() -> None:
    descriptor = PedigreeDescriptor(
        family_id="fabricated-collapse",
        members=tuple(_member(i) for i in range(1, 6)),
        relationships=(_edge(1, 3), _edge(1, 4), _edge(3, 5), _edge(4, 5)),
    )
    assert validate_pedigree(descriptor).family_id == "fabricated-collapse"


def _genotype(alt_count: int) -> tuple[str, str]:
    return (("A", "A"), ("A", "G"), ("G", "G"))[alt_count]


@pytest.mark.parametrize(("first", "second", "child"), tuple(product(range(3), repeat=3)))
def test_exhaustive_trio_truth_table(first: int, second: int, child: int) -> None:
    result = evaluate_site(
        variant_id="fabricated-variant",
        child_member_id="fabricated-child",
        child_genotype=_genotype(child),
        parent_genotypes=(_genotype(first), _genotype(second)),
        relationship_assertion_ids=("edge-1", "edge-2"),
        genotype_record_ids=("gt-1", "gt-2", "gt-3"),
        reference="A",
        alternate="G",
        chromosome="1",
    )
    possible = {
        a + b
        for a in ({0} if first == 0 else {1} if first == 2 else {0, 1})
        for b in ({0} if second == 0 else {1} if second == 2 else {0, 1})
    }
    expected = (
        CompatibilityStatus.CONSISTENT if child in possible else CompatibilityStatus.INCONSISTENT
    )
    assert result.compatibility == expected
    assert len(result.assignments) == len(
        set((x.parent_1_allele, x.parent_2_allele) for x in result.assignments)
    )


def test_duo_unknown_parent_is_not_assumed_reference() -> None:
    result = evaluate_site(
        variant_id="fabricated-variant",
        child_member_id="fabricated-child",
        child_genotype=("G", "G"),
        parent_genotypes=(("A", "A"),),
        relationship_assertion_ids=("edge-1",),
        genotype_record_ids=("gt-1", "gt-2"),
        reference="A",
        alternate="G",
        chromosome="1",
    )
    assert result.compatibility == CompatibilityStatus.INCONSISTENT
    assert result.transmission_status == TransmissionStatus.NO_COMPATIBLE_TRANSMISSION


def _evaluate_special(
    *,
    child_genotype: tuple[str, str] | None = ("A", "G"),
    chromosome: str = "1",
    conflicting: bool = False,
    unsupported_ploidy: bool = False,
) -> CompatibilityStatus:
    return evaluate_site(
        variant_id="fabricated-variant",
        child_member_id="fabricated-child",
        child_genotype=child_genotype,
        parent_genotypes=(("A", "A"), ("G", "G")),
        relationship_assertion_ids=("edge-1", "edge-2"),
        genotype_record_ids=("gt-1",),
        reference="A",
        alternate="G",
        chromosome=chromosome,
        conflicting=conflicting,
        unsupported_ploidy=unsupported_ploidy,
    ).compatibility


def test_missing_conflicting_and_unsupported_are_explicit() -> None:
    assert _evaluate_special(child_genotype=None) == CompatibilityStatus.INDETERMINATE_MISSING
    assert _evaluate_special(conflicting=True) == CompatibilityStatus.INDETERMINATE_CONFLICTING
    assert _evaluate_special(chromosome="X") == CompatibilityStatus.UNSUPPORTED_LOCUS
    assert _evaluate_special(unsupported_ploidy=True) == CompatibilityStatus.UNSUPPORTED_PLOIDY
    assert (
        _evaluate_special(child_genotype=("A", "T"))
        == CompatibilityStatus.UNSUPPORTED_GENOTYPE_REPRESENTATION
    )
