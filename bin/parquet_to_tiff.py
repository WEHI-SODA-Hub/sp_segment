#!/usr/bin/env python
'''
Module      : parquet_to_tiff
Description : Converts a Cellpose Parquet file containing segmentation
              geometries into a label-image TIFF file.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import typer
import numpy as np
import pyarrow.parquet as pq
import rasterio
from rasterio.features import rasterize
from tifffile import TiffFile, imwrite
from shapely import wkb
from typing import Annotated


def get_image_dimensions(tiff_path: Path) -> tuple[int, int]:
    """
    Extracts image dimensions from the OME-TIFF metadata.
    """
    with TiffFile(tiff_path) as tiff:
        ome_metadata = tiff.ome_metadata

    root = ET.fromstring(ome_metadata)
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

    # Extract dimensions from Pixels element
    pixel_data = root.findall('.//ome:Pixels', ns)[0]
    size_x = int(pixel_data.attrib['SizeX'])
    size_y = int(pixel_data.attrib['SizeY'])

    return size_x, size_y


def main(
    parquet: Annotated[Path, typer.Argument(
        help="Path to Cellpose Parquet input file."
    )],
    tiff_path: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file that the segmentation was run on."
    )],
    geometry_colname: Annotated[str, typer.Option(
        help="Column name containing WKB geometries."
    )] = 'geometry',
):
    """
    Converts a Cellpose Parquet file containing segmentation geometries
    into a label-image TIFF file.
    """
    table = pq.read_table(parquet)
    df = table.to_pandas()

    # Load geometries
    geometries = df[geometry_colname].apply(
        lambda geom: wkb.loads(geom)
    ).values

    # Create incrementing IDs for each geometry (1-based as 0 is background)
    ids = [id for id in range(1, len(geometries) + 1)]

    # Get X and Y dimensions from the TIFF file
    (size_x, size_y) = get_image_dimensions(tiff_path)

    # Rasterize: (geometry, value)
    shapes = ((geom, fid) for geom, fid in zip(geometries, ids))
    bounds = (
        0, 0, size_x, size_y
    )
    mask = rasterize(
        shapes=shapes,
        out_shape=(size_y, size_x),  # Note: rasterio uses (height, width)
        fill=0,
        transform=rasterio.transform.from_bounds(
            *bounds,
            width=size_x,
            height=size_y
        ),
        dtype=np.uint16,
        all_touched=True  # Fill shapes
    )

    # Write TIFF output
    imwrite(sys.stdout.buffer, np.flipud(mask))


if __name__ == "__main__":
    typer.run(main)
