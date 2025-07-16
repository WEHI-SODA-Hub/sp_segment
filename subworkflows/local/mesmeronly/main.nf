include { MESMERSEGMENT as MESMERWC  } from '../../../modules/local/mesmersegment/main.nf'
include { MESMERSEGMENT as MESMERNUC } from '../../../modules/local/mesmersegment/main.nf'
include { CELLMEASUREMENT            } from '../../../modules/local/cellmeasurement/main.nf'

workflow MESMERONLY {

    take:
    ch_mesmeronly // channel: segmentation parameters

    main:

    ch_versions = Channel.empty()

    ch_mesmeronly.map {
        sample,
        run_backsub,
        run_mesmer,
        run_cellpose,
        tiff,
        nuclear_channel,
        membrane_channels,
        mesmer_combine_method,
        mesmer_level,
        mesmer_maxima_threshold,
        mesmer_interior_threshold,
        mesmer_maxima_smooth,
        mesmer_min_nuclei_area,
        mesmer_remove_border_cells,
        mesmer_pixel_expansion,
        mesmer_padding,
        skip_measurements -> [
            sample,
            tiff,
            nuclear_channel,
            membrane_channels,
            mesmer_combine_method,
            mesmer_level,
            mesmer_maxima_threshold,
            mesmer_interior_threshold,
            mesmer_maxima_smooth,
            mesmer_min_nuclei_area,
            mesmer_remove_border_cells,
            mesmer_pixel_expansion,
            mesmer_padding
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
    ch_mesmeronly
        .join(MESMERNUC.out.segmentation_mask)
        .join(MESMERWC.out.segmentation_mask)
        .map {
            sample,
            run_backsub,
            run_mesmer,
            run_cellpose,
            tiff,
            nuclear_channel,
            membrane_channels,
            mesmer_combine_method,
            mesmer_level,
            mesmer_maxima_threshold,
            mesmer_interior_threshold,
            mesmer_maxima_smooth,
            mesmer_min_nuclei_area,
            mesmer_remove_border_cells,
            mesmer_pixel_expansion,
            mesmer_padding,
            skip_measurements,
            nuclear_mask,
            whole_cell_mask -> [
                sample,
                tiff,
                nuclear_mask,
                whole_cell_mask,
                skip_measurements
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
