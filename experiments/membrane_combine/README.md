# Membrane channel combination comparison

Compares raw per-pixel maximum (the pipeline's method) with rescaling each
membrane channel before taking the maximum or mean, and with a single marker,
across Cellpose, Mesmer and CellSAM.

The method, results and conclusion are in
[docs/background/membrane-channels.md](../../docs/background/membrane-channels.md).

- `make_inputs.py`: builds the membrane inputs and the samplesheet from a
  config.
- `analyse.py`: scores the pipeline's masks and draws the overlay figures.
- `config_brain_run5.json`: test 1, one COMET brain crop.
- `config_csl.json`: test 2, five crops from two COMET atherosclerosis images,
  with manually labelled cells as ground truth.
- `brain_scores.csv`, `brain_contrast.csv`: `analyse.py` output for test 1
  (2026-10-14).
- `brain_per_marker.csv`: per-marker measures for test 1 from an earlier
  version of `analyse.py`. Its nuclei-inside-a-cell values differ from
  `brain_scores.csv` by up to 0.7 points. That version used blob sigmas of
  5–15 pixels; the current one uses 5.05–15.15 pixels (radius 2–6 µm at
  0.28 µm/pixel). Both find 153 nuclei.
- `csl_scores.csv`, `csl_contrast.csv`, `csl_gt_by_type.csv`: `analyse.py`
  output for test 2 (2026-10-14).

Usage:

```bash
python make_inputs.py --config config_csl.json --outdir <inputs>
# run the pipeline on <inputs>/samplesheet.csv with --deepcell_cache_dir set
python analyse.py --config config_csl.json --inputs <inputs> --results <outdir> --out <dir>
```

Both scripts import from `bin/`, so run them from a clone of the repository.
They need numpy, tifffile, zarr, scipy, scikit-image and matplotlib.
