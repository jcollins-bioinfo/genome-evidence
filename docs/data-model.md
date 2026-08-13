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

## M1 source records

`RawGenotypeObservation` is an immutable pre-normalization record containing the unchanged marker, chromosome and genotype tokens, positive source position, source line, pseudonymous sample, run ID, call state, and explicitly lexical categories. It intentionally has no REF, ALT, canonical genotype, annotation, dosage, or interpretation fields. `SourceMetadata` retains exact-file hash/size, comments, explicit versus overridden build provenance, counts, parser/package version, timestamp, and run. `AssayQCSummary` and `QCFinding` provide aggregate and record-referenced descriptive QC.
