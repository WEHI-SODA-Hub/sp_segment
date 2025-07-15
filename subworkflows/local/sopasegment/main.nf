/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
include { SOPACELLPOSE as SOPACELLPOSENUCLEAR     } from '../sopacellpose/main.nf'
include { SOPACELLPOSE as SOPACELLPOSEWHOLECELL   } from '../sopacellpose/main.nf'
include { PARQUETTOTIFF as PARQUETTOTIFFNUCLEAR   } from '../../../modules/local/parquettotiff/main.nf'
include { PARQUETTOTIFF as PARQUETTOTIFFWHOLECELL } from '../../../modules/local/parquettotiff/main.nf'
include { CELLMEASUREMENT                         } from '../../../modules/local/cellmeasurement/main.nf'

process SOPACONVERT {
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7'
        : 'docker.io/quentinblampey/sopa:2.0.7'}"

    input:
    tuple val(meta), path(tiff)

    output:
    tuple val(meta), path("*.zarr"), emit: spatial_data

    script:
    def args = task.ext.args ?: ''
    """
    sopa convert \\
        ${args} \\
        --sdata-path ${meta.id}.zarr \\
        --technology ${params.technology} \\
        ${tiff}
    """
}

process SOPAPATCHIFYIMAGE {
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7'
        : 'docker.io/quentinblampey/sopa:2.0.7'}"

    input:
    tuple val(meta), path(zarr)

    output:
    tuple val(meta), path("*.zarr/.sopa_cache/patches_file_image"), path("*.zarr/shapes/image_patches"), emit: patches

    script:
    def args = task.ext.args ?: ''
    """
    sopa patchify image \\
        ${args} \\
        ${zarr} \\
        --patch-width-pixel ${params.patch_width_pixel} \\
        --patch-overlap-pixel ${params.patch_overlap_pixel}
    """
}

workflow SOPASEGMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = Channel.empty()

    ch_sopa.map { meta, tiff, nuclear_channel, membrane_channels, skip_measurements ->
        [ meta, tiff ]
    }.set { ch_convert }

    //
    // Run SOPA convert to convert tiff to zarr format
    //
    SOPACONVERT(
        ch_convert
    )

    //
    // Run SOPA patchify to create image patches
    //
    SOPAPATCHIFYIMAGE(
        SOPACONVERT.out.spatial_data
    )

    // Create a channel for each patch
    SOPAPATCHIFYIMAGE.out.patches
        .join( SOPACONVERT.out.spatial_data, by: 0 )
        .map { meta, patches_file_image, image_patches, zarr ->
            [ meta, zarr, patches_file_image.text.trim().toInteger() ] }
        .flatMap { meta, zarr, n_patches ->
            (0..<n_patches).collect { index -> [ meta, zarr, index, n_patches ] } }
        .combine(ch_sopa, by: 0)
        .map { meta, zarr, index, n_patches, tiff, nuclear_channel, membrane_channels, skip_measurements ->
            [ meta, zarr, index, n_patches, nuclear_channel.first(), membrane_channels.first() ]
        }.set { ch_cellpose }

    //
    // Run SOPA with cellpose for nuclear segmentation
    //
    SOPACELLPOSENUCLEAR(
        ch_cellpose.map { meta, zarr, index, n_patches, nuclear_channel, membrane_channels ->
            [ meta, zarr, index, n_patches, nuclear_channel, [] ]
        }, // remove membrane channels for nuclear segmentation
        SOPACONVERT.out.spatial_data
    )

    //
    // Run SOPA with cellpose for whole-cell segmentation
    //
    SOPACELLPOSEWHOLECELL(
        ch_cellpose,
        SOPACONVERT.out.spatial_data
    )

    //
    // Convert nuclear segmentation parquet to tiff
    //
    PARQUETTOTIFFNUCLEAR(
        SOPACELLPOSENUCLEAR.out.boundaries
            .join(ch_sopa, by: 0)
            .map { meta, boundaries, tiff, nuc_chan, mem_chans, skip_measure ->
                [ meta, boundaries, tiff, 'nuclear' ]
            }
    )

    //
    // Convert whole-cell segmentation parquet to tiff
    //
    PARQUETTOTIFFWHOLECELL(
        SOPACELLPOSEWHOLECELL.out.boundaries
            .join(ch_sopa, by: 0)
            .map { meta, boundaries, tiff, nuc_chan, mem_chans, skip_measure ->
                [ meta, boundaries, tiff, 'whole-cell' ]
            }
    )

    //
    // Create a channel for cell measurement
    //
    PARQUETTOTIFFNUCLEAR.out.tiff
        .join(PARQUETTOTIFFWHOLECELL.out.tiff, by: 0)
        .join(ch_sopa, by: 0)
        .map {
            meta,
            nuclear_tiff,
            wholecell_tiff,
            tiff,
            nuc_chan,
            mem_chans,
            skip_measurements -> [
                meta,
                tiff,
                nuclear_tiff,
                wholecell_tiff,
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
    zarr                 = SOPACONVERT.out.spatial_data         // channel: [ val(meta), *.zarr ]
    nuclear_boundaries   = SOPACELLPOSENUCLEAR.out.boundaries   // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]
    wholecell_boundaries = SOPACELLPOSEWHOLECELL.out.boundaries // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                      // channel: [ versions.yml ]
}
