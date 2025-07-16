include { BACKGROUNDSUBTRACT } from '../backgroundsubtract/main.nf'
include { MESMER_SEGMENT     } from '../mesmer_segment/main.nf'

workflow MESMER_SEGMENT_WBACKSUB {

    take:
    ch_mesmer_wbacksub

    main:

    ch_versions = Channel.empty()

    //
    // Run background subtraction
    //
    BACKGROUNDSUBTRACT(
        ch_mesmer_wbacksub.map {
            sample,
            _run_backsub,
            _run_mesmer,
            _run_cellpose,
            tiff,
            _nuclear_channel,
            _membrane_channels,
            _mesmer_combine_method,
            _mesmer_level,
            _mesmer_maxima_threshold,
            _mesmer_interior_threshold,
            _mesmer_maxima_smooth,
            _mesmer_min_nuclei_area,
            _mesmer_remove_border_cells,
            _mesmer_pixel_expansion,
            _mesmer_padding -> [
                sample,
                tiff
            ]
        }
    )

    // Replace tiff with backsub_tif
    ch_mesmer_wbacksub
        .join( BACKGROUNDSUBTRACT.out.backsub_tif )
        .map {
            sample,
            run_backsub,
            run_mesmer,
            run_cellpose,
            _tiff,
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
            backsub_tiff -> [
                sample,
                run_backsub,
                run_mesmer,
                run_cellpose,
                backsub_tiff,
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

    MESMER_SEGMENT(
        ch_mesmer
    )

    emit:
    annotations      = MESMER_SEGMENT.out.annotations       // channel: [ val(meta), *.geojson ]
    whole_cell_tif   = MESMER_SEGMENT.out.whole_cell_tif    // channel: [ val(meta), *.tiff ]
    nuclear_tif      = MESMER_SEGMENT.out.nuclear_tif       // channel: [ val(meta), *.tiff ]

    versions = ch_versions                                  // channel: [ versions.yml ]
}
