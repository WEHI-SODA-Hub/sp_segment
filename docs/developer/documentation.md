# Documentation

The documentation site is built with
[MkDocs](https://www.mkdocs.org/) and
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) from the
Markdown files in `docs/`. The page order is set in `nav` in `mkdocs.yml`.

## Build locally

```bash
python3 -m venv .venv-docs
.venv-docs/bin/pip install -r docs/requirements.txt
.venv-docs/bin/mkdocs serve        # preview at http://127.0.0.1:8000
.venv-docs/bin/mkdocs build --strict
```

`--strict` fails the build on any warning, including broken links. CI uses it.

## Generated pages

The tables in `docs/reference/parameters.md` and
`docs/reference/samplesheet.md` are generated from `nextflow_schema.json` and
`assets/schema_input.json`. Only the text between the
`<!-- BEGIN GENERATED -->` and `<!-- END GENERATED -->` markers is replaced.

After changing either schema, run:

```bash
python3 docs/scripts/generate_reference.py
```

CI runs it with `--check` and fails if the committed pages are out of date.
Prettier is told to skip these two files (`.prettierignore`), so it does not
reformat the generated tables.

## Writing style

- Use short sentences and literal language.
- Give exact parameter names and default values.
- Describe what the pipeline does now. Put history in `CHANGELOG.md`.
- Put implementation detail in the Background or Developer pages, not in the
  how-to and reference pages.
- Use GitHub-style callouts (`> [!NOTE]`, `> [!WARNING]`). They render on
  GitHub and on the site.
- Do not put code blocks inside numbered lists. Prettier and MkDocs indent list
  content differently.
- Link between pages with relative `.md` links, so links work on GitHub and on
  the site. Link to files outside `docs/` with a full GitHub URL.

## Deployment

`.github/workflows/docs.yml` builds the site on every pull request and deploys
it with [mike](https://github.com/jimporter/mike) to the `gh-pages` branch:

| Event             | Deployed as                               |
| ----------------- | ----------------------------------------- |
| push to `dev`     | version `dev`                             |
| release published | version `<tag>`, alias `latest` (default) |

The site is served at <https://wehi-soda-hub.github.io/sp_segment/>.

> [!NOTE]
> A repository admin must enable GitHub Pages once: **Settings > Pages >
> Build and deployment > Source: Deploy from a branch**, branch `gh-pages`,
> folder `/ (root)`. The branch is created by the first deployment.
