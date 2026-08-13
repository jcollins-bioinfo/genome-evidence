# Canonical private workspace

The provider-neutral workspace defaults in personal Colab to `/content/drive/MyDrive/genome-evidence-private`; `GENOME_EVIDENCE_WORKSPACE` overrides it. Manual uploads belong only in `inputs/inbox/23andme/`. Initialization is idempotent and rejects Git checkouts. Exact bytes are copied to `inputs/raw/23andme/<sha256>/genome.txt`; cache is evictable while inputs, references, registries, runs, reports, and exports are durable.

Use only pseudonyms such as `subject-0001`. The JSON registries contain workspace-relative paths; no symlink is used for latest. Compute in ephemeral local storage, copy into a unique Drive destination, verify every hash there, then write `COMPLETED.json` last. Incomplete runs remain under `_incomplete` and are not registered. Drive FUSE rename atomicity is never assumed.

Colab executes on a Google-managed VM. Users who require no cloud processing must run the same path locally. Personal mode never substitutes synthetic data, and production marker, FASTA, evidence, M5, and M6 bundles must be installed and validated explicitly.
