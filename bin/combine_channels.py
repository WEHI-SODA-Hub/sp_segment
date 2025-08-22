#!/usr/bin/env python
'''
Module      : combine_channels
Description : Takes an OME-TIFF containing N channels and returns a 2 channel tiff
              containing a nuclear channel and membrane channel that is componsed
              of one or more channels from the input tiff, usign either the product
              or max of the intensities.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import re
import sys
import xml.etree.ElementTree as ET
from enum import Enum
from pathlib import Path
from typing import Annotated, List

import typer
import numpy as np
from tifffile import TiffFile, imwrite
from xarray import DataArray, concat


class CombineMethod(str, Enum):
    PROD = "prod"
    MAX = "max"


def get_pixels_tag(xml_str: str) -> ET.Element:
    """
    Parses the OME-XML string and returns the Pixels tag.
    """
    root = ET.fromstring(xml_str)
    ns = {'ome': root.tag.split('}')[0].strip('{')}
    image_tag = root.find('ome:Image', ns)

    if image_tag is None:
        raise ValueError("No Image tag found in the XML.")

    pixels_tag = image_tag.find('ome:Pixels', ns)
    if pixels_tag is None:
        raise ValueError("No Pixels tag found in the Image tag.")

    return pixels_tag


def extract_channel_names(xml_str) -> List[str]:
    """
    Extracts all channel 'Name' attributes from OME-XML in the order they
    appear. Returns a list of channel names by index.
    """
    pixels_tag = get_pixels_tag(xml_str)

    root = ET.fromstring(xml_str)
    ns = {'ome': root.tag.split('}')[0].strip('{')}
    channel_tags = pixels_tag.findall('ome:Channel', ns)

    channel_names = [ch.get('Name', 'Unknown') for ch in channel_tags]
    return channel_names


def tiff_to_xarray(tiffPath: Path) -> DataArray:
    """
    Takes a TIFF and converts it to an xarray with relevant axis,
    coordinate and metadata attached. Supports MIBI TIFF and OME-TIFF.
    """
    channel_names: list[str] = []
    attrs: dict[str, float] = {}
    #: List of channels, each of which are 2D
    channels = []

    with TiffFile(tiffPath) as tiff:
        first_page = tiff.pages[0]

        # OME-TIFF: channel info in first page only
        channel_names = extract_channel_names(first_page.description)
        for page in tiff.pages:
            channels.append(page.asarray())

    return DataArray(data=channels, dims=["C", "Y", "X"],
                     coords={"C": channel_names}, attrs=attrs)


def combine_channels(array: DataArray, channels: List[str], combined_name: str,
                     combine_method: CombineMethod) -> DataArray:
    """
    Combines multiple channels into a single channel using the specified method
    (prod or max). Adds the combined channel to the array.
    """

    if len(channels) == 1:
        return array

    combined = array.sel(C=channels)

    if combine_method == CombineMethod.MAX:
        combined = combined.max(dim="C")
    elif combine_method == CombineMethod.PROD:
        combined = combined.prod(dim="C")

    combined = combined.expand_dims("C").assign_coords(C=[combined_name])

    return concat([array, combined], dim="C")


def update_ome_xml(original_xml, width, height, channels, channel_names):
    """
    Update existing OME-XML metadata with new channel information.
    """
    if original_xml is None or not original_xml.startswith('<?xml'):
        raise ValueError("Invalid OME-XML metadata provided.")

    try:
        original_xml = re.sub(r'<\?xml[^>]+\?>', '<?xml version="1.0"?>', original_xml)
        root = ET.fromstring(original_xml)

        namespace = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
        pixels = root.find('.//ome:Pixels', namespace)

        if pixels is not None:
            # Update dimensions
            pixels.set('SizeX', str(width))
            pixels.set('SizeY', str(height))
            pixels.set('SizeC', str(channels))

            # Update/set PlaneCount to match channel number
            tiff_data_elements = pixels.findall('ome:TiffData', namespace)
            if tiff_data_elements:
                for tiff_data in tiff_data_elements:
                    tiff_data.set('PlaneCount', str(channels))
            else:
                tiff_data = ET.SubElement(pixels, 'TiffData')
                tiff_data.set('PlaneCount', str(channels))

            for channel in pixels.findall('ome:Channel', namespace):
                pixels.remove(channel)

            for c in range(channels):
                channel_id = f"Channel:{c}"
                channel_attrs = {
                    'ID': channel_id,
                    'Name': channel_names[c],
                    'SamplesPerPixel': '1'
                }
                ET.SubElement(pixels, 'Channel', channel_attrs)

        xml_str = '<?xml version="1.0" encoding="UTF-8"?>\n' + \
            ET.tostring(root, encoding='utf-8').decode('utf-8')

        # Fix ns0 namespace prefixes if they exists
        xml_str = re.sub(r'<ns0:', '<', xml_str)
        xml_str = re.sub(r'</ns0:', '</', xml_str)
        xml_str = re.sub(r'xmlns:ns0', 'xmlns', xml_str)

        return xml_str

    except ET.ParseError as e:
        raise ValueError(f"Failed to parse OME-XML metadata: {e}")


def main(
    tiff: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file."
    )],
    nuclear_channel: Annotated[str, typer.Option(
        help="Name of the nuclear channel."
    )],
    membrane_channel: Annotated[List[str], typer.Option(
        help="Name(s) of the membrane channels (can be repeated)"
             "Ensure that channels with spaces are quoted.")
    ],
    combine_method: Annotated[CombineMethod, typer.Option(
        help="Method to use for combining channels (prod or max).")
    ] = CombineMethod.PROD,
):
    full_array = tiff_to_xarray(tiff)

    # Combine channels and prepare image array
    combined_membrane_channel = "combined_membrane" \
        if len(membrane_channel) > 1 else membrane_channel[0]
    full_array = combine_channels(full_array, membrane_channel,
                                  combined_membrane_channel,
                                  CombineMethod(combine_method))

    # Extract the nuclear and membrane channels
    output_array = full_array.sel(
        C=[nuclear_channel, combined_membrane_channel]
    ).values.astype(np.uint16)

    with TiffFile(tiff) as tif:
        ome_metadata = tif.pages[0].description

    # Update OME-XML metadata
    c, height, width = output_array.shape
    updated_metadata = update_ome_xml(ome_metadata, width, height, c,
                                      [nuclear_channel, combined_membrane_channel])

    imwrite(sys.stdout.buffer, output_array,
            photometric='minisblack',
            metadata={'axes': 'CYX'},
            description=updated_metadata.encode('utf-8'))


if __name__ == "__main__":
    typer.run(main)
