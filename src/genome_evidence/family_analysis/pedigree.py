"""Strict pedigree parsing and directed ancestry-graph validation."""

from pathlib import Path

from pydantic import ValidationError

from .models import PedigreeDescriptor


def load_pedigree(path: Path) -> PedigreeDescriptor:
    """Load and structurally validate a local JSON pedigree descriptor.

    Parameters
    ----------
    path:
        Checked local JSON file. Its member ordering never affects the result.

    Returns
    -------
    PedigreeDescriptor
        Canonically ordered, structurally valid declared assertions.
    """
    try:
        descriptor = PedigreeDescriptor.model_validate_json(path.read_bytes())
    except (OSError, ValidationError) as error:
        raise ValueError(
            "pedigree descriptor is missing, invalid, or uses an unknown schema"
        ) from error
    return validate_pedigree(descriptor)


def validate_pedigree(descriptor: PedigreeDescriptor) -> PedigreeDescriptor:
    """Reject ambiguous bindings and directed parent cycles without inferring kinship."""
    if not descriptor.members or not descriptor.relationships:
        raise ValueError("pedigree must contain members and an analyzable parent relationship")
    member_ids = [item.member_id for item in descriptor.members]
    subject_ids = [item.subject_id for item in descriptor.members]
    if len(member_ids) != len(set(member_ids)):
        raise ValueError("duplicate pedigree member ID")
    if len(subject_ids) != len(set(subject_ids)):
        raise ValueError("one subject cannot be bound as multiple pedigree members")
    run_bindings = [(str(item.m2_run), item.subject_id) for item in descriptor.members]
    if len({path for path, _ in run_bindings}) != len(run_bindings):
        raise ValueError("one M2 run cannot be bound to multiple subjects")
    known = set(member_ids)
    edge_pairs: set[tuple[str, str]] = set()
    assertion_ids: set[str] = set()
    parents: dict[str, set[str]] = {}
    children: dict[str, set[str]] = {member: set() for member in known}
    for edge in descriptor.relationships:
        if edge.assertion_id in assertion_ids:
            raise ValueError("duplicate relationship assertion ID")
        assertion_ids.add(edge.assertion_id)
        pair = (edge.parent_member_id, edge.child_member_id)
        if pair in edge_pairs:
            raise ValueError("duplicate biological-parent edge")
        edge_pairs.add(pair)
        if edge.parent_member_id not in known or edge.child_member_id not in known:
            raise ValueError("relationship references an absent member")
        if edge.parent_member_id == edge.child_member_id:
            raise ValueError("self-parent relationship is invalid")
        parents.setdefault(edge.child_member_id, set()).add(edge.parent_member_id)
        children[edge.parent_member_id].add(edge.child_member_id)
    if any(len(items) > 2 for items in parents.values()):
        raise ValueError("initial M9 supports at most two declared biological parents")

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(member: str) -> None:
        if member in visiting:
            raise ValueError("directed biological-parent ancestry cycle")
        if member in visited:
            return
        visiting.add(member)
        for child in sorted(children[member]):
            visit(child)
        visiting.remove(member)
        visited.add(member)

    for member in sorted(known):
        visit(member)
    return PedigreeDescriptor(
        family_id=descriptor.family_id,
        members=tuple(sorted(descriptor.members, key=lambda item: item.member_id)),
        relationships=tuple(
            sorted(
                descriptor.relationships,
                key=lambda item: (item.child_member_id, item.parent_member_id, item.assertion_id),
            )
        ),
    )
