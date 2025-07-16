/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap        } from 'plugin/nf-schema'
include { BACKGROUNDSUBTRACT      } from '../subworkflows/local/backgroundsubtract'
include { MESMER_SEGMENT_WBACKSUB } from '../subworkflows/local/mesmer_segment_wbacksub'
include { MESMER_SEGMENT          } from '../subworkflows/local/mesmer_segment'
include { SOPA_SEGMENT            } from '../subworkflows/local/sopa_segment'
include { SOPA_SEGMENT_WBACKSUB   } from '../subworkflows/local/sopa_segment_wbacksub'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_spatialproteomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPATIALPROTEOMICS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    //
    // Construct channel for background subtraction/segmentation workflow
    //
    ch_samplesheet.map {
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
        mesmer_padding -> [
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
            mesmer_padding
        ]
    }.branch { it ->
        backsub_only: it[1].contains(true) && !it[2].contains(true) &&
                        !it[3].contains(true) // run_backsub true, run_mesmer false, run_cellpose false
        backsub_mesmer: it[1].contains(true) && it[2].contains(true) // run_backsub true, run_mesmer true
        mesmer_only: !it[1].contains(true) && it[2].contains(true) // run_backsub false, run_mesmer true
    }.set { ch_segmentation_samplesheet }

    // TODO: setup test data compatible with background subtraction
    //
    // Run the BACKGROUNDSUBTRACT subworkflow for samples that ONLY require
    // background subtraction (no segmentation)
    //
    BACKGROUNDSUBTRACT(
        ch_segmentation_samplesheet.backsub_only.map {
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

    //
    // Run the MESMER_SEGMENT_WBACKSUB subworkflow for samples that require
    // background subtraction and mesmer segmentation
    //
    MESMER_SEGMENT_WBACKSUB(
        ch_segmentation_samplesheet.backsub_mesmer
    )

    //
    // Run MESMER_SEGMENT subworkflow for samples that ONLY require mesmer segmentation
    //
    MESMER_SEGMENT(
        ch_segmentation_samplesheet.mesmer_only
    )

    //
    // Construct channel for only CELLPOSE subworkflow
    //
    ch_samplesheet.filter {
        it[3].contains(true) // run_cellpose true for sample
    }.map {
        sample,
        run_backsub,
        _run_mesmer,
        _run_cellpose,
        tiff,
        nuclear_channel,
        membrane_channels,
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
            run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels
        ]
    }.branch { it ->
        with_backsub: it[1].contains(true)// run_backsub true
        no_backsub: !it[1].contains(true) // run_backsub false
    }.set { ch_cellpose_samplesheet }

    //
    // Run CELLPOSE subworkflow for samples that require background subtraction
    //
    SOPA_SEGMENT_WBACKSUB(
        ch_cellpose_samplesheet.with_backsub.map { sample,
            _run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels ->
            [ sample, tiff, nuclear_channel, membrane_channels ]
        }
    )

    //
    // Run CELLPOSE subworkflow for samples that ONLY require cellpose segmentation
    //
    SOPA_SEGMENT(
        ch_cellpose_samplesheet.no_backsub.map { sample,
            _run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels ->
            [ sample, tiff, nuclear_channel, membrane_channels ]
        }
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  ''  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_collated_versions     // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
