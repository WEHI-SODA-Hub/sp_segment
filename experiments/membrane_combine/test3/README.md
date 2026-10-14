# Test 3: objects without a nucleus, missed cells and outlines

Results are described in
[docs/background/membrane-channels.md](../../../docs/background/membrane-channels.md#test-3-objects-without-a-nucleus-missed-cells-and-outlines).

The scripts read the test 2 inputs (`make_inputs.py` with `config_csl.json`),
the test 2 pipeline results, and the production cells of the five crops
exported from the CSL QuPath project. The paths are WEHI scratch paths set at
the top of each script. They were run in the pipeline's sopa-Cellpose
container (2026-10-15).

Run in this order:

1. `score_artefacts.py <outdir>`: rasterises the production cells of each crop
   to `<crop>_production_mask.npy`, which the other scripts read.
2. `low_dapi.py`: low-DAPI cells per input → `low_dapi.csv`.
3. `missed_cells.py`: raw-max cells with clear DAPI that the other inputs lack
   → `missed_vs_raw_max.csv`.
4. `segment_variants.py`, `segment_variants2.py` (GPU): Cellpose setting
   variants on each crop (diameter, cell probability and flow thresholds, sopa
   preprocessing, a 3-channel input).
5. `score_shape.py`: reference outlines from each labelled cell's marker and
   nucleus, and IoU, recall and precision for every mask →
   `shape_summary.csv`, `shape_by_type.csv`.
6. `overall.py`: cells, low-DAPI cells, ground truth and shape scores per mask
   → `overall.csv`.

In the tables, `(pipeline)` masks come from the test 2 pipeline run. The other
masks come from Cellpose run directly on each crop with the same settings
except the one named; `raw_max cellpose base` is that run with default
settings.
