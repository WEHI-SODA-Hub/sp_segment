# Run several methods on one image

Each samplesheet row runs at most one segmentation method, and each `sample`
value must be unique. To run two or three methods on one image, add one row per
method with a different sample name and the same `tiff`:

```csv
sample,run_backsub,run_mesmer,run_cellpose,run_cellsam,tiff,nuclear_channel,membrane_channels
slide1_cellpose,true,false,true,false,/data/slide1.ome.tiff,DAPI,CD45:Pan-Cytokeratin
slide1_mesmer,true,true,false,false,/data/slide1.ome.tiff,DAPI,CD45:Pan-Cytokeratin
slide1_cellsam,true,false,false,true,/data/slide1.ome.tiff,DAPI,CD45:Pan-Cytokeratin
```

Each row writes its own files, for example
`cellmeasurement/slide1_cellpose.geojson.gz` and
`cellmeasurement/slide1_mesmer.geojson.gz`.

Background subtraction runs once per row, so this example runs it three times.
To run it once, run a first pass with `run_backsub` only, then use
`backsub/<sample>.ome.tif` as the `tiff` for the segmentation rows with
`run_backsub = false`.

The pipeline stops at startup if two rows share a `sample` value, or if a row
sets more than one of `run_mesmer`, `run_cellpose` and `run_cellsam`.
