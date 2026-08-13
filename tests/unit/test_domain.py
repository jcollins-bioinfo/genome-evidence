from datetime import UTC, datetime

import pytest
from pydantic import ValidationError

from genome_evidence.domain.enums import CallStatus, PhasingStatus
from genome_evidence.domain.inference import GenotypeInference
from genome_evidence.domain.observations import GenotypeObservation
from genome_evidence.domain.samples import Sample
from genome_evidence.domain.subjects import Subject
from genome_evidence.domain.variants import Variant


def sample() -> Sample:
    return Sample(
        sample_id="synthetic-sample",
        subject=Subject(subject_id="synthetic-subject"),
        source="synthetic",
        assay_type="array",
        genome_build="GRCh37",
        source_identifier="fixture",
        ingested_at=datetime(2026, 1, 1, tzinfo=UTC),
        source_checksum="sha256:synthetic",
    )


def test_models_instantiate_and_types_remain_distinct() -> None:
    variant = Variant(assembly="GRCh38", chromosome="1", position=101, reference="A", alternate="G")
    observation = GenotypeObservation(
        observation_id="obs-1",
        sample=sample(),
        source_marker_id="marker-1",
        source_build="GRCh37",
        source_chromosome="1",
        source_position=100,
        original_genotype="AG",
        original_alleles=("A", "G"),
        call_status=CallStatus.CALLED,
        observation_method="synthetic array",
    )
    inference = GenotypeInference(
        inference_id="inf-1",
        sample_id=sample().sample_id,
        variant=variant,
        genotype_probabilities={"A/A": 0.1, "A/G": 0.8, "G/G": 0.1},
        dosage=1.0,
        imputation_quality=0.9,
        reference_panel="synthetic-panel",
        reference_panel_version="1",
        phasing_status=PhasingStatus.UNPHASED,
        run_id="run-1",
    )
    assert not isinstance(observation, GenotypeInference)
    assert not isinstance(inference, GenotypeObservation)


def test_invalid_coordinate_fails_validation() -> None:
    with pytest.raises(ValidationError):
        Variant(assembly="GRCh38", chromosome="1", position=0, reference="A", alternate="G")


def test_missing_is_not_homozygous_reference() -> None:
    common = dict(
        sample=sample(),
        source_build="GRCh37",
        source_chromosome="3",
        source_position=303,
        observation_method="synthetic array",
    )
    missing = GenotypeObservation(
        observation_id="missing",
        source_marker_id="missing",
        original_genotype=None,
        original_alleles=(),
        call_status=CallStatus.NO_CALL,
        **common,
    )
    reference = GenotypeObservation(
        observation_id="reference",
        source_marker_id="reference",
        original_genotype="CC",
        original_alleles=("C", "C"),
        call_status=CallStatus.CALLED,
        **common,
    )
    assert missing.call_status != reference.call_status
    assert missing.original_alleles != reference.original_alleles


def test_source_observation_is_immutable() -> None:
    observation = GenotypeObservation(
        observation_id="obs-immutable",
        sample=sample(),
        source_marker_id="marker",
        source_build="GRCh37",
        source_chromosome="1",
        source_position=1,
        original_genotype=None,
        original_alleles=(),
        call_status=CallStatus.NO_CALL,
        observation_method="synthetic array",
    )
    with pytest.raises(ValidationError):
        observation.original_genotype = "AA"  # type: ignore[misc]
