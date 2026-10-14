# Membrane channel combination comparison

Compares raw per-pixel maximum (the pipeline's method) with rescaling each
membrane channel before taking the maximum or mean, and with a single marker,
across Cellpose, Mesmer and CellSAM.

The method, results and conclusion are in
[docs/background/membrane-channels.md](../../docs/background/membrane-channels.md).

- `make_inputs.py`: builds the membrane inputs and the samplesheet.
- `analyse.py`: scores the pipeline's masks and draws the overlay figure.
- `results_run5_crop.csv`, `contrast_run5_crop.csv`: results for the test image
  (COMET brain crop, 2026-10-14).

Both scripts import from `bin/`, so run them from a clone of the repository.
They need numpy, tifffile, scipy, scikit-image and matplotlib.
