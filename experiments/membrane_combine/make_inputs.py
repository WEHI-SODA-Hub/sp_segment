#!/usr/bin/env python3
"""
Build the membrane-channel images for the membrane combination comparison.

Writes one 2-channel OME-TIFF per combination method, each holding the nuclear
channel and one membrane channel, plus a samplesheet that runs every image
through Cellpose, Mesmer and CellSAM. Because each image already has a single
membrane channel, the pipeline's COMBINECHANNELS step passes it through
unchanged, and the segmenters see exactly the image built here.

Methods:
  raw_max    per-pixel maximum of the raw channels (the pipeline's method)
  norm_max   each channel rescaled (median -> 0, 99.5th percentile -> 1,
             negative values clipped to 0), then the per-pixel maximum
  norm_mean  the same rescaling, then the per-pixel mean
  single     one reference marker on its own

Usage:
  python make_inputs.py --tiff <background-subtracted.ome.tif> --outdir <dir>
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))
from combine_channels import build_ome_xml, source_physical_size, to_uint16  # noqa: E402
from tiff_channels import get_channel_names, read_channels  # noqa: E402

MARKERS = ["VCW_CD45", "VCW_MHCII", "VCW_Iba1_AIF1", "VCW_GFAP_diluted", "VCW_Integrin_b1"]
SINGLE = "VCW_Integrin_b1"
NUCLEAR = "DAPI"
SEGMENTERS = ["cellpose", "mesmer", "cellsam"]


def rescale(plane):
    lo, hi = np.percentile(plane, [50, 99.5])
    return np.clip((plane.astype(np.float64) - lo) / (hi - lo), 0, None)


def to_full_range(plane):
    """Linear map to 0..65535 (every segmenter rescales its input anyway)."""
    peak = float(plane.max())
    return to_uint16(plane * (65535.0 / peak)) if peak > 0 else to_uint16(plane)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiff", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    names = get_channel_names(args.tiff)
    nuclear = read_channels(args.tiff, [names.index(NUCLEAR)])[0]
    planes = read_channels(args.tiff, [names.index(m) for m in MARKERS])
    normalised = np.stack([rescale(p) for p in planes])

    membranes = {
        "raw_max": to_uint16(planes.max(axis=0)),
        "norm_max": to_full_range(normalised.max(axis=0)),
        "norm_mean": to_full_range(normalised.mean(axis=0)),
        "single": to_uint16(planes[MARKERS.index(SINGLE)]),
    }

    physical_size = source_physical_size(Path(args.tiff))
    rows = []
    for method, membrane in membranes.items():
        path = out / f"{method}.ome.tif"
        data = np.stack([to_uint16(nuclear), membrane])
        xml = build_ome_xml(data.shape[2], data.shape[1], [NUCLEAR, "membrane"], physical_size)
        tifffile.imwrite(path, data, photometric="minisblack", metadata=None,
                         description=xml.encode("utf-8"))
        for seg in SEGMENTERS:
            rows.append({
                "sample": f"{method}_{seg}", "run_backsub": "false",
                "run_mesmer": str(seg == "mesmer").lower(),
                "run_cellpose": str(seg == "cellpose").lower(),
                "run_cellsam": str(seg == "cellsam").lower(),
                "tiff": str(path.resolve()), "nuclear_channel": NUCLEAR,
                "membrane_channels": "membrane",
            })
        print(f"wrote {path}")

    with open(out / "samplesheet.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out / 'samplesheet.csv'} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
