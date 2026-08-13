# Initial data model

| Entity | Purpose | Key relationships |
|---|---|---|
| `Subject` | Project-local pseudonymous person | Has samples; no real name required |
| `Sample` | Assay/source instance and checksum | Belongs to a subject; produces observations |
| `GenotypeObservation` | Immutable, original source call or explicit missing state | Belongs to sample; may have mappings |
| `Variant` | Assembly/chromosome/position/REF/ALT canonical allele | Target of mappings and inferences; rsID is optional metadata, not identity |
| `ObservationVariantMapping` | Derived source-to-canonical mapping | Links observation ID to optional variant; records strand, liftover, ambiguity, confidence, and run |
| `GenotypeInference` | Probabilistic genotype/dosage | Links sample and variant; includes panel, quality, phase, and run; never substitutes for an observation |
| `EvidenceAssertion` | Versioned directed relationship between typed entities | Carries evidence category, interpretation status, source/accession/version/retrieval date, and run |
| `RunProvenance` | Reproducibility envelope | Identifies inputs, hashes, software/git/configuration, references, transformation, and times |

Missingness is encoded through observation call status with absent alleles. It is structurally different from a called observation whose alleles equal the reference. Assertion endpoints are generic typed references so relationships such as variant→phenotype, variant→gene, gene→disease, variant→polygenic score, and drug→pharmacogene fit a relational representation.

DuckDB is the initial analytical source of truth, backed where practical by immutable/versioned Parquet artifacts and source checksums. Database state should be reconstructible from source and derived artifacts. Pydantic validation protects API boundaries; database constraints and schemas will mature with ingestion milestones.

## M2 derived records

`ObservationMapping` links a deterministic M1 observation reference and normalization run to an optional canonical target, explicit status/reason, strand action, liftover state/candidates, allele correspondence, and REF validation. `CanonicalGenotype` exists only for defensible called alleles and preserves haploid/diploid ploidy. `ResourceIdentity` records local resource provenance. These types bridge rather than populate the older M0 `GenotypeObservation`, whose fields M1 cannot justify.

## M1 source records

`RawGenotypeObservation` is an immutable pre-normalization record containing the unchanged marker, chromosome and genotype tokens, positive source position, source line, pseudonymous sample, run ID, call state, and explicitly lexical categories. It intentionally has no REF, ALT, canonical genotype, annotation, dosage, or interpretation fields. `SourceMetadata` retains exact-file hash/size, comments, explicit versus overridden build provenance, counts, parser/package version, timestamp, and run. `AssayQCSummary` and `QCFinding` provide aggregate and record-referenced descriptive QC.

## M3 external records

`ExternalSourceSnapshot` identifies exact release bytes and XML metadata. `ExternalVariantRepresentation` retains source identifiers and alleles. `ExternalAssertion` retains a submitted SCV or aggregate VCV claim, source vocabulary, logical accession/version, content fingerprint, and snapshot instance identity. `AssertionRelationship` keeps submission-to-aggregate relationships explicit. `VariantEvidenceLink` connects a representation—not a genotype or assertion—to an M2 `variant_id`, with explicit matched/unmatched/ambiguous/incompatible/unsupported outcomes.

## M4 review-routing records

`VariantEvidenceProfile` joins exact-linked source evidence to a canonical variant while preserving genotype state, SCV/VCV assertions, classification dimensions, terms, conditions, and unresolved assessments. `ClinicalReviewCandidate` carries a non-clinical review band and explicit eligibility. `CandidateAssertionLink`, `PriorityRationale`, and `PrioritizationExclusion` normalize traceability; `PolicyIdentity` captures exact policy bytes and canonical parsed configuration.

## M5 model outputs

`MarkerAlignment` retains every model marker and its exact used, missing, excluded, or discordant state. `ProjectionCoordinate`, `ReferenceGroupDistance`, `ReferenceNeighbor`, and `ReferenceSupportEvaluation` are bundle-dependent model outputs. `ProjectionSensitivityReplicate` records deterministic leave-one-chromosome-out sensitivity.

## M6 inference records

M6 keeps target alignments/exclusions, observed statistical phase inference, imputed genotype inference, official variant quality, masked internal-consistency rows, and chromosome state/checkpoint rows as distinct typed artifacts. Every scientific identifier is content/lineage derived. Official fields retain ALT orientation and header-declared Number/Type semantics. Zero-row tables still have explicit schemas.
