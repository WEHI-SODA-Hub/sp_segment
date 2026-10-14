#!/usr/bin/env python3
"""
Build the membrane-channel images for the membrane combination comparison.

For each region in the config, writes:
  <region>_source.ome.tif    the nuclear channel and the membrane markers, used
                             by analyse.py
  <region>_<method>.ome.tif  the nuclear channel and one membrane channel built
                             with <method>
and one samplesheet that runs every <region>_<method> image through each
segmenter. Each image already has a single membrane channel, so the pipeline's
COMBINECHANNELS step passes it through unchanged and the segmenters see exactly
the image built here.

Methods:
  raw_max    per-pixel maximum of the raw channels (the pipeline's method)
  norm_max   each channel rescaled (median -> 0, 99.5th percentile -> 1,
             values below 0 set to 0), then the per-pixel maximum
  norm_mean  the same rescaling, then the per-pixel mean
  single     one marker on its own: the config's "single" marker, or with
             "auto" the marker that wins the most pixels in raw_max

Config (JSON):
  {
    "nuclear": "DAPI",
    "markers": ["CD45", ...],
    "single": "auto",
    "segmenters": ["cellpose", "mesmer", "cellsam"],
    "regions": [
      {"name": "crop1", "tiff": "/path/image.ome.tif", "y": 0, "x": 0, "size": 2048},
      {"name": "whole", "tiff": "/path/small.ome.tif"}
    ],
    "ground_truth": "/path/labels.csv"      (optional, used by analyse.py)
  }
A region without y/x/size uses the whole image.

Usage:
  python make_inputs.py --config <config.json> --outdir <dir>
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from combine_channels import build_ome_xml, source_physical_size, to_uint16  # noqa: E402
from tiff_channels import get_channel_names, read_channels  # noqa: E402

METHODS = ["raw_max", "norm_max", "norm_mean", "single"]


def rescale(plane):
    lo, hi = np.percentile(plane, [50, 99.5])
    return np.clip((plane.astype(np.float64) - lo) / (hi - lo), 0, None)


def to_full_range(plane):
    """Linear map to 0..65535 (every segmenter rescales its input anyway)."""
    peak = float(plane.max())
    return to_uint16(plane * (65535.0 / peak)) if peak > 0 else to_uint16(plane)


def read_region(path, indices, region):
    """(C, H, W) planes for the channel indices, cropped when the region says so."""
    if "size" not in region:
        return read_channels(path, indices)
    import zarr

    y, x, size = region["y"], region["x"], region["size"]
    with tifffile.imread(path, aszarr=True) as store:
        z = zarr.open(store, mode="r")
        if hasattr(z, "keys") and not hasattr(z, "shape"):
            z = z[sorted(z.keys())[0]]  # pyramidal: level 0 is full resolution
        return np.stack([np.asarray(z[i, y:y + size, x:x + size]) for i in indices])


def write_ome(path, planes, names, physical_size):
    xml = build_ome_xml(planes.shape[2], planes.shape[1], names, physical_size)
    tifffile.imwrite(path, planes, photometric="minisblack", metadata=None,
                     description=xml.encode("utf-8"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    nuclear_name, markers = cfg["nuclear"], cfg["markers"]
    segmenters = cfg.get("segmenters", ["cellpose", "mesmer", "cellsam"])
    rows, singles = [], {}
    for region in cfg["regions"]:
        name = region["name"]
        names = get_channel_names(region["tiff"])
        indices = [names.index(nuclear_name)] + [names.index(m) for m in markers]
        stack = read_region(region["tiff"], indices, region)
        nuclear, planes = to_uint16(stack[0]), stack[1:]
        physical_size = source_physical_size(Path(region["tiff"]))
        write_ome(out / f"{name}_source.ome.tif", np.stack([to_uint16(p) for p in stack]),
                  [nuclear_name] + markers, physical_size)

        raw_max = planes.max(axis=0)
        single = cfg.get("single", "auto")
        if single == "auto":
            wins = np.bincount(np.argmax(planes, axis=0).ravel(), minlength=len(markers))
            single = markers[int(np.argmax(wins))]
        singles[name] = single
        normalised = np.stack([rescale(p) for p in planes])
        membranes = {
            "raw_max": to_uint16(raw_max),
            "norm_max": to_full_range(normalised.max(axis=0)),
            "norm_mean": to_full_range(normalised.mean(axis=0)),
            "single": to_uint16(planes[markers.index(single)]),
        }
        for method in METHODS:
            path = out / f"{name}_{method}.ome.tif"
            write_ome(path, np.stack([nuclear, membranes[method]]), [nuclear_name, "membrane"],
                      physical_size)
            for seg in segmenters:
                rows.append({
                    "sample": f"{name}_{method}_{seg}", "run_backsub": "false",
                    "run_mesmer": str(seg == "mesmer").lower(),
                    "run_cellpose": str(seg == "cellpose").lower(),
                    "run_cellsam": str(seg == "cellsam").lower(),
                    "tiff": str(path.resolve()), "nuclear_channel": nuclear_name,
                    "membrane_channels": "membrane",
                })
        print(f"{name}: wrote source and {len(METHODS)} inputs (single = {single})")

    with open(out / "samplesheet.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (out / "singles.json").write_text(json.dumps(singles, indent=2))
    print(f"wrote {out / 'samplesheet.csv'} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
