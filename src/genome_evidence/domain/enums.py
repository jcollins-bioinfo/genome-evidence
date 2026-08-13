"""Closed vocabularies for domain boundaries; not clinical grading systems."""

from enum import StrEnum


class EvidenceType(StrEnum):
    DIRECT_OBSERVATION = "direct_observation"
    IMPUTED_GENOTYPE = "imputed_genotype"
    DERIVED_HAPLOTYPE = "derived_haplotype"
    EXTERNAL_ASSOCIATION = "external_association"
    EXTERNAL_CLINICAL_ASSERTION = "external_clinical_assertion"
    MEASURED_PHENOTYPE = "measured_phenotype"
    SELF_REPORTED_PHENOTYPE = "self_reported_phenotype"
    FAMILY_HISTORY = "family_history"
    MODEL_OUTPUT = "model_output"


class InterpretationStatus(StrEnum):
    ESTABLISHED = "established"
    SUPPORTED = "supported"
    PROVISIONAL = "provisional"
    SPECULATIVE = "speculative"
    INDETERMINATE = "indeterminate"


class CallStatus(StrEnum):
    CALLED = "called"
    NO_CALL = "no_call"
    UNASSAYED = "unassayed"


class MappingStatus(StrEnum):
    MAPPED = "mapped"
    UNMAPPED = "unmapped"
    AMBIGUOUS = "ambiguous"
    FAILED = "failed"


class StrandTransformation(StrEnum):
    NONE = "none"
    COMPLEMENT = "complement"
    REVERSE_COMPLEMENT = "reverse_complement"
    UNKNOWN = "unknown"


class LiftoverStatus(StrEnum):
    NOT_REQUIRED = "not_required"
    SUCCESS = "success"
    FAILED = "failed"
    AMBIGUOUS = "ambiguous"


class MappingConfidence(StrEnum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    UNKNOWN = "unknown"


class PhasingStatus(StrEnum):
    UNPHASED = "unphased"
    PHASED = "phased"
    PARTIALLY_PHASED = "partially_phased"
