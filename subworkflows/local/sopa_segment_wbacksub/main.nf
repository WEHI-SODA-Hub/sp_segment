include { BACKGROUNDSUBTRACT } from '../backgroundsubtract/main.nf'
include { SOPA_SEGMENT       } from '../sopa_segment/main.nf'

workflow SOPA_SEGMENT_WBACKSUB {

    take:
    ch_sopa_wbacksub

    main:

    ch_versions = Channel.empty()

    //
    // Run the BACKGROUNDSUBTRACT subworkflow for samples that ONLY require
    // background subtraction (no segmentation)
    //
    BACKGROUNDSUBTRACT(
        ch_sopa_wbacksub.map {
            sample,
            tiff,
            _nuclear_channel,
            _membrane_channels -> [
                sample,
                tiff
            ]
        }
    )

    // Replace tiff with backsub_tif
    ch_sopa_wbacksub
        .join( BACKGROUNDSUBTRACT.out.backsub_tif )
        .map { sample,
            _tiff,
            nuclear_channel,
            membrane_channels,
            backsub_tiff ->
            [ sample, backsub_tiff, nuclear_channel, membrane_channels ]
        }.set { ch_sopa }

    SOPA_SEGMENT(
        ch_sopa
    )

    emit:
    annotations = SOPA_SEGMENT.out.annotations // channel: [ val(meta), *.geojson ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
