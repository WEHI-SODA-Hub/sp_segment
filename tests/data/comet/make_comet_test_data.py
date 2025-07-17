#!/usr/bin/env python3
'''
Create a synthetic multichannel image with OME-XML metadata to mimic a COMET
TIFF file.

The image is licensed under MIT and was obtained from nf-core/test-datasets:
https://github.com/nf-core/test-datasets/raw/refs/heads/modules/data/imaging/segmentation/nuclear_image.tif

Code was mostly written by Claude Sonnet 4.
'''
import xml.etree.ElementTree as ET
from tifffile import TiffFile, TiffWriter
import numpy as np


def create_synthetic_multichannel_image(
        input_path='nuclear_image.tif',
        output_path='synthetic_multichannel.tif'
):
    with TiffFile(input_path) as tiff:
        original_data = tiff.asarray()
        original_metadata = tiff.ome_metadata

    root = ET.fromstring(original_metadata)
    pixels = root.find(
            './/{http://www.openmicroscopy.org/Schemas/OME/2016-06}Pixels'
    )

    size_x = int(pixels.attrib['SizeX'])
    size_y = int(pixels.attrib['SizeY'])

    print(f"Original image dimensions: {size_x} x {size_y}")

    # Create synthetic 3-channel data
    # Channel 0: DAPI (original nuclear image)
    channel_0 = original_data.astype(np.uint16)

    # Channel 1: TRITC (constant background)
    tritc_background_level = 500  # Constant background value
    channel_1 = np.full((size_y, size_x), tritc_background_level,
                        dtype=np.uint16)

    # Channel 2: CD45 with TRITC background (background + original signal)
    channel_2 = (channel_1 + channel_0).astype(np.uint16)

    # Stack channels together (C, Y, X format)
    multichannel_data = np.stack([channel_0, channel_1, channel_2], axis=2)

    ome_xml = create_ome_metadata(size_x, size_y, num_channels=3)

    with TiffWriter(output_path) as tiff:
        tiff.write(
            multichannel_data,
            metadata={'axes': 'YXC'},
            description=ome_xml,
            # Add TIFF tags to indicate this is an OME-TIFF
            extratags=[
                ('ImageDescription', 's', 0, ome_xml, True),  # OME-XML in ImageDescription
            ]
        )
    print(f"Synthetic multichannel image saved to: {output_path}")
    return multichannel_data, ome_xml


def create_ome_metadata(size_x, size_y, num_channels=3):
    """Create OME-XML metadata with Plane and ChannelPriv elements"""

    # Channel configurations
    channels_config = [
        {'name': 'DAPI', 'color': '-16776961', 'exposure_time': '25.0', 'led_current': '20.0', 'background': None},
        {'name': 'TRITC', 'color': '-65536', 'exposure_time': '250.0', 'led_current': '1700.0', 'background': None},
        {'name': 'CD45', 'color': '-256', 'exposure_time': '250.0', 'led_current': '1700.0', 'background': 'TRITC'}
    ]

    # Build OME-XML structure
    ome_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
     Creator="synthetic_image_generator"
     xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">
    <Image ID="Image:0" Name="SyntheticMultichannel">
        <Pixels ID="Pixels:0:0" DimensionOrder="XYZCT" SizeC="{num_channels}" SizeT="1" SizeX="{size_x}" SizeY="{size_y}" SizeZ="1" Type="uint16">'''

    ome_xml += '''
        <TiffData/>
    '''

    # Add Channel elements
    for i, config in enumerate(channels_config):
        ome_xml += f'''
            <Channel ID="Channel:{i}" Color="{config['color']}" Name="{config['name']}" SamplesPerPixel="1"/>'''

    # Add Plane elements
    for i in range(num_channels):
        ome_xml += f'''
            <Plane ExposureTime="{channels_config[i]['exposure_time']}" ExposureTimeUnit="ms" TheC="{i}" TheT="0" TheZ="0"/>'''

    ome_xml += '''
        </Pixels>
    </Image>
    <StructuredAnnotations>
        <XMLAnnotation ID="Annotation:0">
            <Value>
                <PrivateFields>'''

    # Add ChannelPriv elements
    for i, config in enumerate(channels_config):
        fluorescence_channel = config['background'] if config['background'] else config['name']
        ome_xml += f'''
                    <ChannelPriv ID="Channel:{i}" CycleID="0" LedCurrent="{config['led_current']}" LedCurrentUnit="mA" SensorGain="0.0" FluorescenceChannel="{fluorescence_channel}"/>'''

    ome_xml += '''
                </PrivateFields>
            </Value>
        </XMLAnnotation>
    </StructuredAnnotations>
</OME>'''

    return ome_xml


if __name__ == "__main__":
    # Generate synthetic multichannel image
    data, metadata = create_synthetic_multichannel_image()

    print(f"\nGenerated image shape: {data.shape}")
    print("Channels: DAPI (0), TRITC (1), CD45 with TRITC background (2)")
