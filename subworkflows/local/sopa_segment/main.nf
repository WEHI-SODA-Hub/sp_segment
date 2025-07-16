/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
include { SOPA_CELLPOSE as SOPA_CELLPOSENUCLEAR   } from '../sopa_cellpose/main.nf'
include { SOPA_CELLPOSE as SOPA_CELLPOSEWHOLECELL } from '../sopa_cellpose/main.nf'
include { SOPA_CONVERT                            } from '../../../modules/local/sopa/convert/main.nf'
include { SOPA_PATCHIFYIMAGE                      } from '../../../modules/local/sopa/patchifyimage/main.nf'
include { PARQUETTOTIFF as PARQUETTOTIFFNUCLEAR   } from '../../../modules/local/parquettotiff/main.nf'
include { PARQUETTOTIFF as PARQUETTOTIFFWHOLECELL } from '../../../modules/local/parquettotiff/main.nf'
include { CELLMEASUREMENT                         } from '../../../modules/local/cellmeasurement/main.nf'

workflow SOPA_SEGMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = Channel.empty()

    ch_sopa.map {
        meta,
        tiff,
        _nuclear_channel,
        _membrane_channels -> [ meta, tiff ]
    }.set { ch_convert }

    //
    // Run SOPA convert to convert tiff to zarr format
    //
    SOPA_CONVERT(
        ch_convert
    )

    //
    // Run SOPA patchify to create image patches
    //
    SOPA_PATCHIFYIMAGE(
        SOPA_CONVERT.out.spatial_data
    )

    // Create a channel for each patch
    SOPA_PATCHIFYIMAGE.out.patches
        .join( SOPA_CONVERT.out.spatial_data, by: 0 )
        .map { meta, patches_file_image, _image_patches, zarr ->
            [ meta, zarr, patches_file_image.text.trim().toInteger() ] }
        .flatMap { meta, zarr, n_patches ->
            (0..<n_patches).collect { index -> [ meta, zarr, index, n_patches ] } }
        .combine(ch_sopa, by: 0)
        .map { meta, zarr, index, n_patches, _tiff, nuclear_channel, membrane_channels ->
            [ meta, zarr, index, n_patches, nuclear_channel.first(), membrane_channels.first() ]
        }.set { ch_cellpose }

    //
    // Run SOPA with cellpose for nuclear segmentation
    //
    SOPA_CELLPOSENUCLEAR(
        ch_cellpose.map { meta, zarr, index, n_patches, nuclear_channel, _membrane_channels ->
            [ meta, zarr, index, n_patches, nuclear_channel, [] ]
        }, // remove membrane channels for nuclear segmentation
        SOPA_CONVERT.out.spatial_data
    )

    //
    // Run SOPA with cellpose for whole-cell segmentation
    //
    SOPA_CELLPOSEWHOLECELL(
        ch_cellpose,
        SOPA_CONVERT.out.spatial_data
    )

    //
    // Convert nuclear segmentation parquet to tiff
    //
    PARQUETTOTIFFNUCLEAR(
        SOPA_CELLPOSENUCLEAR.out.boundaries
            .join(ch_sopa, by: 0)
            .map { meta, boundaries, tiff, _nuc_chan, _mem_chans ->
                [ meta, boundaries, tiff, 'nuclear' ]
            }
    )

    //
    // Convert whole-cell segmentation parquet to tiff
    //
    PARQUETTOTIFFWHOLECELL(
        SOPA_CELLPOSEWHOLECELL.out.boundaries
            .join(ch_sopa, by: 0)
            .map { meta, boundaries, tiff, _nuc_chan, _mem_chans ->
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
            _nuc_chan,
            _mem_chans -> [
                meta,
                tiff,
                nuclear_tiff,
                wholecell_tiff
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
    zarr                 = SOPA_CONVERT.out.spatial_data         // channel: [ val(meta), *.zarr ]
    nuclear_boundaries   = SOPA_CELLPOSENUCLEAR.out.boundaries   // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]
    wholecell_boundaries = SOPA_CELLPOSEWHOLECELL.out.boundaries // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                      // channel: [ versions.yml ]
}
