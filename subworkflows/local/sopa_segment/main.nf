include { SOPA_SEGMENT_COMPARTMENT as SOPA_SEGMENT_NUCLEAR   } from '../sopa_segment_compartment/main.nf'
include { SOPA_SEGMENT_COMPARTMENT as SOPA_SEGMENT_WHOLECELL } from '../sopa_segment_compartment/main.nf'
include { CELLMEASUREMENT                                    } from '../../../modules/local/cellmeasurement/main.nf'

workflow SOPA_SEGMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = Channel.empty()

    //
    // Run segmentation for nuclear compartment
    //
    SOPA_SEGMENT_NUCLEAR(
        ch_sopa.map {
            meta,
            tiff,
            nuclear_channel,
            _membrane_channels -> [
                meta,
                tiff,
                nuclear_channel,
                [""] // no membrane channels for nuclear segmentation
            ]
        },
        'nuclear'
    )

    //
    // Run segmentation for whole-cell compartment
    //
    SOPA_SEGMENT_WHOLECELL(
        ch_sopa,
        'whole-cell'
    )

    //
    // Create a channel for cell measurement
    //
    SOPA_SEGMENT_NUCLEAR.out.tiff
        .join(SOPA_SEGMENT_WHOLECELL.out.tiff, by: 0)
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
    annotations = CELLMEASUREMENT.out.annotations   // channel: [ val(meta), *.geojson ]

    versions = ch_versions                          // channel: [ versions.yml ]
}
