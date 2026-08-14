# Agent guide

- Inspect existing work and preserve unrelated changes.
- Read `docs/epistemic-contract.md` before changing domain semantics; update architecture docs when semantics change.
- Never commit real genomic or clinical data. Tests use clearly synthetic fixtures only, and logging must not dump genotypes.
- Maintain provenance and the hard observed-versus-inferred distinction. Do not add speculative biological logic without explicit scope.
- Add or update tests for behavioral changes.
- Before completion run `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src`, and `uv run pytest`.
- Report failures and verification gaps; never conceal them.
- For every change, assess package-version and changelog impact, public API/schema
  compatibility, root README impact, and documentation/type/test impact.
- Update the root README whenever installation, CLI, notebooks, milestones,
  architecture, outputs, dependencies, limitations, privacy, or workflow changes.
