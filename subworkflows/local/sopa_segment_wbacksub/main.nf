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
            nuclear_channel,
            membrane_channels -> [
                sample,
                tiff
            ]
        }
    )

    // Replace tiff with backsub_tif
    ch_sopa_wbacksub
        .join( BACKGROUNDSUBTRACT.out.backsub_tif )
        .map { sample,
            tiff,
            nuclear_channel,
            membrane_channels,
            backsub_tiff ->
            [ sample, backsub_tiff, nuclear_channel, membrane_channels ]
        }.set { ch_sopa }

    SOPA_SEGMENT(
        ch_sopa
    )

    emit:
    zarr                 = SOPA_SEGMENT.out.zarr                 // channel: [ val(meta), *.zarr ]
    nuclear_boundaries   = SOPA_SEGMENT.out.nuclear_boundaries   // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]
    wholecell_boundaries = SOPA_SEGMENT.out.wholecell_boundaries // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                       // channel: [ versions.yml ]
}
