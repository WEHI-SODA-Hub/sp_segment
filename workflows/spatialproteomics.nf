/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap       } from 'plugin/nf-schema'
include { BACKGROUNDSUBTRACT     } from '../subworkflows/local/backgroundsubtract'
include { BACKSUBMESMER          } from '../subworkflows/local/backsubmesmer'
include { MESMERONLY             } from '../subworkflows/local/mesmeronly'
include { SOPA_SEGMENT           } from '../subworkflows/local/sopa_segment'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spatialproteomics_pipeline'

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
        mesmer_padding,
        skip_measurements -> [
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
            skip_measurements
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
            ]
        }
    )

    //
    // Run the BACKSUBMESMER subworkflow for samples that require
    // background subtraction and mesmer segmentation
    //
    BACKSUBMESMER(
        ch_segmentation_samplesheet.backsub_mesmer
    )

    //
    // Run MESMERONLY subworkflow for samples that ONLY require mesmer segmentation
    //
    MESMERONLY(
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
            skip_measurements
        ]
    }.set { ch_cellpose_samplesheet }

    // TODO: add sopa segment with bg subtract workflow
    //
    // Run CELLPOSE subworkflow for samples that ONLY require cellpose segmentation
    //
    SOPA_SEGMENT(
        ch_cellpose_samplesheet
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
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
