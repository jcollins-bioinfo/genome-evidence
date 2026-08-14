# Versioning and releases

Milestone, package, artifact schema, external source release, and task/prompt
revision are independent identities. `pyproject.toml` is the sole package-version
declaration; runtime and CLI reporting use installed distribution metadata.

Before 1.0, a backwards-compatible milestone feature increments `0.y.0`, a
backwards-compatible correction increments `0.y.z`, and incompatible removals need
a documented deprecation/migration path. Version 1.0 requires intentionally stable
public APIs, schemas, CLI contracts, migrations, and compatibility commitments.
Schema readers accept documented schema IDs, never infer compatibility from the
package version, and fail closed on unknown majors.

## Authorized release checklist

1. Classify API/schema compatibility and select a PEP 440 version.
2. Update `pyproject.toml`, changelog, README, migrations, and compatibility docs.
3. Run all locked quality, notebook, build, and clean-wheel checks.
4. Confirm package/runtime/CLI/wheel version agreement and inspect artifacts.
5. Obtain explicit release authorization; review credentials, signing, and supply chain.
6. Create the signed tag, GitHub release, and publication only in that authorized task.

M8 does not automate or perform tagging, GitHub releases, signing, or publication.
