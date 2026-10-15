"""
Tests for bin/extract_markers.py: the exposure time and background channel of
each marker must stay on that marker, because BACKSUB scales and subtracts
with them.
"""

from extract_markers import ome_channel_meta

NS = 'xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"'


def _by_name(meta):
    return {cm["marker_name"]: cm for cm in meta}


def test_background_follows_channel_id_not_position():
    # ChannelPriv elements listed in the reverse order of the channels.
    xml = f"""<OME {NS}><Image ID="Image:0"><Pixels ID="Pixels:0">
        <Channel ID="Channel:0" Name="DAPI"/>
        <Channel ID="Channel:1" Name="TRITC"/>
        <Channel ID="Channel:2" Name="CD45"/>
        <Plane TheC="0" ExposureTime="25"/>
        <Plane TheC="1" ExposureTime="250"/>
        <Plane TheC="2" ExposureTime="300"/>
      </Pixels></Image>
      <StructuredAnnotations><XMLAnnotation><Value>
        <ChannelPriv ID="Channel:2" FluorescenceChannel="TRITC"/>
        <ChannelPriv ID="Channel:1" FluorescenceChannel="TRITC"/>
        <ChannelPriv ID="Channel:0" FluorescenceChannel="DAPI"/>
      </Value></XMLAnnotation></StructuredAnnotations></OME>"""
    meta = _by_name(ome_channel_meta(xml))
    assert meta["CD45"]["background"] == "TRITC"
    assert meta["TRITC"]["background"] == ""
    assert meta["DAPI"]["background"] == ""
    assert [meta[n]["exposure"] for n in ("DAPI", "TRITC", "CD45")] == [25.0, 250.0, 300.0]


def test_one_channel_per_image_keeps_each_exposure():
    # Each Image holds one channel, so every Plane has TheC=0.
    images = "".join(
        f'<Image ID="Image:{i}"><Pixels ID="Pixels:{i}">'
        f'<Channel ID="Channel:{i}:0" Name="{name}"/>'
        f'<Plane TheC="0" ExposureTime="{exp}"/></Pixels></Image>'
        for i, (name, exp) in enumerate([("DAPI", 25), ("CD3", 100), ("CD8", 200)])
    )
    meta = ome_channel_meta(f"<OME {NS}>{images}</OME>")
    assert [(cm["marker_name"], cm["exposure"]) for cm in meta] == [
        ("DAPI", 25.0), ("CD3", 100.0), ("CD8", 200.0)
    ]


def test_other_ome_namespace():
    xml = """<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2015-01">
      <Image ID="Image:0"><Pixels ID="Pixels:0">
        <Channel ID="Channel:0" Name="DAPI"/><Plane TheC="0" ExposureTime="25"/>
      </Pixels></Image></OME>"""
    assert ome_channel_meta(xml) == [
        {"index": "Channel:0", "marker_name": "DAPI", "exposure": 25.0}
    ]
