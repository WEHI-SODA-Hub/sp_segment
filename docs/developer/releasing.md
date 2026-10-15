# Releasing

Releases are made from `dev` to `main`. Tags have no `v` prefix (for example
`0.5.0`).

1. On a branch from `dev`, set the new version:
   `nf-core pipelines bump-version <version>`. This updates
   `manifest.version` in `nextflow.config` and the other version strings.
2. In `CHANGELOG.md`, rename `## Unreleased` to `## v<version> - <YYYY-MM-DD>`
   and add a new empty `## Unreleased` above it.
3. Open a pull request to `dev`, then a pull request from `dev` to `main`.
   On pull requests to `main`, lint runs with `--release`.
4. Before merging, run the test profiles affected by the release on a GPU node.
   See [Testing](testing.md).
5. After merging, create a GitHub release from `main` with the tag `<version>`
   and the changelog section as the release notes.

Publishing the release:

- runs CI on the release
- deploys the documentation as version `<version>` with the alias `latest`
  (see [Documentation](documentation.md#deployment))
- creates a new version of the Zenodo record, if the Zenodo GitHub integration
  is enabled for the repository

The pipeline DOI cited in the README is
[10.5281/zenodo.17103183](https://doi.org/10.5281/zenodo.17103183).
