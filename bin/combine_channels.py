#!/usr/bin/env python
'''
Module      : combine_channels
Description : Takes a multichannel TIFF (OME-TIFF including OPAL OME/QPTIFF,
              MIBI TIFF, or ImageJ TIFF) and writes a 2-channel OME-TIFF to
              stdout: the nuclear channel, and a membrane channel that is
              either one input channel or several combined with max or
              product.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import sys
import xml.etree.ElementTree as ET
from enum import Enum
from pathlib import Path
from typing import Annotated, List, Optional

import numpy as np
import typer
from tifffile import TiffFile, imwrite

from tiff_channels import (
    get_channel_names,
    ome_physical_size,
    read_channels,
)

UINT16_MAX = np.iinfo(np.uint16).max
COMBINED_MEMBRANE = "combined_membrane"


class CombineMethod(str, Enum):
    PROD = "prod"
    MAX = "max"


def to_uint16(data: np.ndarray) -> np.ndarray:
    '''
    Cast to uint16 without wrapping: values are rounded and clipped to
    0..65535. A plain astype() wraps negative floats (e.g. from background
    subtraction) to large positive values.
    '''
    if data.dtype == np.uint16:
        return data
    if np.issubdtype(data.dtype, np.floating):
        data = np.rint(data)
    return np.clip(data, 0, UINT16_MAX).astype(np.uint16)


def combine(planes: np.ndarray, method: CombineMethod) -> np.ndarray:
    '''
    Combine a (C, Y, X) stack into one uint16 (Y, X) plane.

    max keeps the input scale. prod is computed in float64, which cannot wrap
    the way an integer product does, and is then rescaled linearly so its
    maximum is 65535 when it exceeds the uint16 range.
    '''
    if method == CombineMethod.MAX:
        return to_uint16(np.max(planes, axis=0))

    combined = np.prod(planes.astype(np.float64), axis=0)
    peak = float(combined.max()) if combined.size else 0.0
    if peak > UINT16_MAX:
        combined *= UINT16_MAX / peak
    return to_uint16(combined)


def channel_index(name: str, channel_names: List[str]) -> int:
    '''Index of a channel name; fails with the list of available names.'''
    try:
        return channel_names.index(name)
    except ValueError:
        raise typer.BadParameter(
            f"channel '{name}' is not in the image. Available channels: {channel_names}"
        )


def build_ome_xml(width: int, height: int, channel_names: List[str],
                  physical_size: Optional[tuple]) -> str:
    '''
    OME-XML for the 2-channel output, built from scratch.

    The source XML is not edited in place: it can use any OME schema version,
    and describes N channels, per-plane TiffData and Plane elements that do
    not apply to the output. Only the pixel size is carried over.
    '''
    ome_ns = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    xsi_ns = "http://www.w3.org/2001/XMLSchema-instance"
    schema_loc = f"{ome_ns} {ome_ns}/ome.xsd"
    ET.register_namespace("", ome_ns)
    ET.register_namespace("xsi", xsi_ns)

    # Clark-notation tag prefixes, built without {{{ patterns that trip nf-core linting
    tag_ome = "{" + ome_ns + "}"
    tag_xsi = "{" + xsi_ns + "}"

    root = ET.Element(tag_ome + "OME")
    root.set(tag_xsi + "schemaLocation", schema_loc)

    image = ET.SubElement(root, tag_ome + "Image", {"ID": "Image:0", "Name": "combined"})
    pixels_attrs = {
        "ID": "Pixels:0",
        "DimensionOrder": "XYZCT",
        "Type": "uint16",
        "SizeX": str(width),
        "SizeY": str(height),
        "SizeZ": "1",
        "SizeC": str(len(channel_names)),
        "SizeT": "1",
    }
    if physical_size is not None:
        size_x, size_y, unit = physical_size
        pixels_attrs.update({
            "PhysicalSizeX": str(size_x),
            "PhysicalSizeY": str(size_y),
            "PhysicalSizeXUnit": unit,
            "PhysicalSizeYUnit": unit,
        })
    pixels = ET.SubElement(image, tag_ome + "Pixels", pixels_attrs)
    # Schema order inside Pixels: Channel elements first, then TiffData.
    for c, name in enumerate(channel_names):
        ET.SubElement(pixels, tag_ome + "Channel", {
            "ID": f"Channel:0:{c}",
            "Name": name,
            "SamplesPerPixel": "1",
        })
    for c in range(len(channel_names)):
        ET.SubElement(pixels, tag_ome + "TiffData", {
            "FirstC": str(c), "FirstT": "0", "FirstZ": "0",
            "IFD": str(c), "PlaneCount": "1",
        })

    return '<?xml version="1.0" encoding="UTF-8"?>\n' + \
        ET.tostring(root, encoding='unicode')


def source_physical_size(tiff: Path) -> Optional[tuple]:
    '''Pixel size from the source OME-XML, or None when it has none.'''
    with TiffFile(tiff) as tif:
        xml = tif.ome_metadata
        if not xml and tif.pages:
            desc = tif.pages[0].description
            xml = desc if isinstance(desc, str) and desc.lstrip().startswith("<") else None
    if not xml:
        return None
    try:
        return ome_physical_size(xml)
    except (ET.ParseError, ValueError):
        return None


def main(
    tiff: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file."
    )],
    nuclear_channel: Annotated[str, typer.Option(
        help="Name of the nuclear channel."
    )],
    membrane_channel: Annotated[List[str], typer.Option(
        help="Name of a membrane channel. Repeat for several channels.")
    ],
    combine_method: Annotated[CombineMethod, typer.Option(
        help="Method to use for combining channels (prod or max).")
    ] = CombineMethod.MAX,
    pixel_size_microns: Annotated[Optional[float], typer.Option(
        help="Pixel size in microns. Written to OME-XML PhysicalSizeX/Y when "
             "the source image has no calibration metadata.")
    ] = None,
):
    channel_names = get_channel_names(tiff)
    nuclear_index = channel_index(nuclear_channel, channel_names)
    membrane_indices = [channel_index(name, channel_names) for name in membrane_channel]

    # Only the nuclear and membrane channels are decoded, not the whole image.
    nuclear = to_uint16(read_channels(tiff, [nuclear_index])[0])
    if len(membrane_indices) > 1:
        membrane = combine(read_channels(tiff, membrane_indices), CombineMethod(combine_method))
        membrane_name = COMBINED_MEMBRANE
    else:
        membrane = to_uint16(read_channels(tiff, membrane_indices)[0])
        membrane_name = membrane_channel[0]

    # Report what was combined so it stays visible in the process log
    # (downstream segmenters only see the single combined membrane channel).
    typer.echo(
        f"[combine_channels] nuclear='{nuclear_channel}'; "
        f"membrane={membrane_channel} -> '{membrane_name}' "
        f"(method={CombineMethod(combine_method).value}, "
        f"n={len(membrane_channel)})",
        err=True,
    )

    output = np.stack([nuclear, membrane])
    _, height, width = output.shape

    physical_size = source_physical_size(tiff)
    if physical_size is None and pixel_size_microns is not None:
        physical_size = (pixel_size_microns, pixel_size_microns, "µm")

    ome_xml = build_ome_xml(width, height, [nuclear_channel, membrane_name], physical_size)

    imwrite(sys.stdout.buffer, output,
            photometric='minisblack',
            metadata=None,
            description=ome_xml.encode('utf-8'))


if __name__ == "__main__":
    typer.run(main)
