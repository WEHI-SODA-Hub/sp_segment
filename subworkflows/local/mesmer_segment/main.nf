include { MESMERSEGMENT as MESMERWC  } from '../../../modules/local/mesmersegment/main.nf'
include { MESMERSEGMENT as MESMERNUC } from '../../../modules/local/mesmersegment/main.nf'
include { CELLMEASUREMENT            } from '../../../modules/local/cellmeasurement/main.nf'

workflow MESMER_SEGMENT {

    take:
    ch_mesmer_segment // channel: segmentation parameters

    main:

    ch_versions = Channel.empty()

    ch_mesmer_segment.map {
        sample,
        _run_backsub,
        _run_mesmer,
        _run_cellpose,
        tiff,
        nuclear_channel,
        membrane_channels -> [
            sample,
            tiff,
            nuclear_channel,
            membrane_channels
        ]
    }.set { ch_mesmer }


    //
    // Run MESMERSEGMENT module on the background subtracted tiff
    // for whole-cell segmentation
    //
    MESMERWC(
        ch_mesmer,
        "whole-cell"
    )

    //
    // Run MESMERSEGMENT module as above, but this time for nuclear segmentation
    //
    MESMERNUC(
        ch_mesmer,
        "nuclear"
    )

    // Create channel for CELLMEASUREMENT input adding the segmentation masks
    ch_mesmer_segment
        .join(MESMERNUC.out.segmentation_mask)
        .join(MESMERWC.out.segmentation_mask)
        .map {
            sample,
            _run_backsub,
            _run_mesmer,
            _run_cellpose,
            tiff,
            _nuclear_channel,
            _membrane_channels,
            nuclear_mask,
            whole_cell_mask -> [
                sample,
                tiff,
                nuclear_mask,
                whole_cell_mask
            ]
        }.set { ch_cellmeasurement }

    //
    // Run CELLMEASUREMENT module on the whole-cell and nuclear segmentation masks
    //
    CELLMEASUREMENT(
        ch_cellmeasurement,
        params.pixel_size_microns
    )

    emit:
    annotations      = CELLMEASUREMENT.out.annotations   // channel: [ val(meta), *.geojson ]
    whole_cell_tif   = MESMERWC.out.segmentation_mask    // channel: [ val(meta), *.tiff ]
    nuclear_tif      = MESMERNUC.out.segmentation_mask   // channel: [ val(meta), *.tiff ]

    versions         = ch_versions                       // channel: [ versions.yml ]
}
