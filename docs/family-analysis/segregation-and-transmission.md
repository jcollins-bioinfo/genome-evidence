# M9 declared-pedigree segregation and transmission evidence

M9 accepts a strict local `genome-evidence-pedigree/v1` JSON descriptor containing a
pseudonymous family ID, unique members, explicit subject-to-completed-M2-run bindings,
and directed biological-parent **assertions** with provenance categories. A valid graph
is not proof of relatedness. Names, birth dates, histories, and narrative fields are
rejected. Each M2 completion marker, manifest, declared checksum, GRCh38 assembly,
`m2-snv-1` algorithm identity, canonical variant key, run lineage, and subject binding
must validate. M6 and every other inferred genotype format have no input interface.

For each child with one or two declared parents, the engine aligns the union of exact
assembly/chromosome/position/REF/ALT canonical identities (represented by the stable M2
variant ID). It never joins by rsID or coordinate alone. Missing sites remain missing;
conflicting calls remain conflicting. Autosomal chromosomes 1–22 and biallelic diploid
allele multisets are the only evaluable scope. Sex chromosomes, PAR, mitochondria,
non-diploid calls, multiallelic/complex representations, CNV/SV/repeats/HLA, and other
unsupported records receive explicit fail-closed statuses.

## Enumeration and result axes

A trio enumerates each distinct ordered pair of alleles transmissible by neutral
`parent_1` and `parent_2` edges, retaining pairs whose unordered multiset equals the
observed child call. A duo enumerates the declared parent's alleles while the unknown
contributor ranges over REF and ALT—never REF alone. At least one assignment is
`CONSISTENT`; none is `INCONSISTENT`. These words are conditional compatibility, not a
de novo call, pedigree verification, non-parentage claim, or causal/clinical result.

Compatibility, informativeness, and site-transmission status are orthogonal. Unique and
ambiguous assignments mean constraints at one site only; they are not haplotype blocks,
recombination, or parent of origin outside the declared-edge context. Missing,
conflicting, unsupported-locus, unsupported-ploidy, and unsupported-representation rows
are never put in the consistent/inconsistent denominator.

The summary's `site_relationship_evidence_count` is the number of canonical-site/child
relationship sets emitted. Each compatibility count uses that same all-row denominator;
no rates are currently reported. If a future report adds a compatibility rate, its
denominator must include only `CONSISTENT + INCONSISTENT`, and must also report excluded
status counts. SNP-array absence means unassayed or missing, never absence from the genome.

An inconsistency can reflect genotyping error, sample mislabelling, inaccurate declared
metadata, normalization/build/strand error, mosaicism, structural or multiallelic
complexity, ploidy assumptions, or a biological event. M9 selects none of these.
Relationship verification would require a separately reviewed kinship model; de novo
calling needs validated sequencing, error/quality and parental-coverage models; sex
chromosomes and SV require locus-specific ploidy/representation contracts; long-range
family phase requires validated phase-block/recombination models; clinical segregation
requires phenotype, penetrance, disease-model and expert-review contracts kept separate.

Outputs are private pseudonymous evidence artifacts with exact assertion, genotype,
manifest, package, Git, configuration and algorithm lineage. Default CLI/notebook output
is aggregate. Cloud notebook execution exposes private artifacts to a cloud VM; local
offline execution is preferred and no regulatory compliance is claimed.
