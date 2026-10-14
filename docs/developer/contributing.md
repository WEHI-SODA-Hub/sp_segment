# Contributing

## Workflow

1. Check the [issues](https://github.com/WEHI-SODA-Hub/sp_segment/issues) for
   existing work on your idea. If there is none, open an issue.
2. Fork the repository, or create a branch if you have write access.
3. Branch from `dev`.
4. Make the change, add or update tests, and update the docs (see
   [Documentation](documentation.md)).
5. Add an entry under `## Unreleased` in `CHANGELOG.md`.
6. Open a pull request against `dev`. `main` is for releases only.

Pull requests are reviewed after the CI checks pass.

## Conventions

- Follow the [nf-core guidelines](https://nf-co.re/docs/guidelines) for
  processes, channels and configuration.
- A new parameter needs a default in the `params` block of `nextflow.config`,
  an entry in `nextflow_schema.json` with the same default, a one-line
  `description` and, if needed, a longer `help_text`. Then regenerate the
  parameter reference with `python3 docs/scripts/generate_reference.py`.
- A new samplesheet column needs an entry in `assets/schema_input.json` with a
  `description`, and a regenerated samplesheet reference.
- Give each process a resource label from `conf/base.config` (`process_low`,
  `process_medium`, `process_high`, `process_multi`, `process_single`,
  `process_gpu`) and add tool arguments in `conf/modules.config` as `ext.args`.
- Each process writes a `versions.yml` with the version of every tool it uses.
  Read versions from the tool, not from a hard-coded string.
- Python scripts in `bin/` need unit tests in `tests/python/`.

## Checks to run before a pull request

```bash
pre-commit run --all-files                       # prettier and editorconfig
pytest tests/python -q                           # bin/ script unit tests
python3 docs/scripts/generate_reference.py --check
mkdocs build --strict                            # needs docs/requirements.txt
nextflow run . -profile test,singularity --outdir results
```

See [Testing](testing.md) for nf-test.

## Code of conduct

See
[CODE_OF_CONDUCT.md](https://github.com/WEHI-SODA-Hub/sp_segment/blob/dev/CODE_OF_CONDUCT.md).
