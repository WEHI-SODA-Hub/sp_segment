# View results in QuPath

The cell measurement GeoJSON can be opened in [QuPath](https://qupath.github.io/)
on top of the image it was made from.

## Which image to open

Open the image the cells were measured on:

- with `run_backsub = true`: `backsub/<sample>.ome.tif`
- otherwise: the `tiff` from the samplesheet

The cell coordinates are in pixels of that image.

## Import the cells

1. Open the image in a QuPath project.
2. Drag `cellmeasurement/<sample>.geojson.gz` onto the image window, or use
   **File > Import objects from file**.
3. Open **Measure > Show detection measurements** to see the measurements.

If QuPath does not accept the `.geojson.gz` file, decompress it first
(`gunzip -k <sample>.geojson.gz`), or set `--gzip_geojson false` to write
uncompressed files.

Measurement names are described in
[Outputs](../output.md#measurement-names).

## Large GeoJSON files

Whole-slide GeoJSON files can be several GB even when compressed, and can be too
large to import through the QuPath interface. The separate import_large_geojson
Nextflow pipeline imports large GeoJSON files into a QuPath project in batch.

## Errors on import

- `StackOverflowError`: some cell outlines have too many vertices. Re-run with
  `--smooth_masks true`, and use `--smooth_method shapely` if that is not
  enough. See [Mask smoothing options](../reference/parameters.md#mask-smoothing-options).
- The cells appear flipped relative to the image (Mesmer only): set
  `--mesmer_transpose_mask true`.
