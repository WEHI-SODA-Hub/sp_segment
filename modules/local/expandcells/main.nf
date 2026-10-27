process EXPANDCELLS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_tifffile_numpy_scipy_pruned:e54488103afb8110"

    input:
    tuple val(meta),
        path(nuclear_mask, stageAs: 'nuclear_mask_input.tiff'),
        path(whole_cell_mask, stageAs: 'whole_cell_mask_input.tiff')

    output:
    tuple val(meta), path("*_expanded.tiff"), emit: mask
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Without nuclear segmentation the nuclear input is a placeholder copy of
    // the whole-cell mask, so it is not passed on.
    // Nuclei are segmented in paired mode, and also when expand_only_without_membrane
    // needs them (segmentNuclei() in utils_nfcore_sp_segment_pipeline).
    def segment_nuclei = (!params.use_whole_cell_only || (params.cell_expansion_microns && params.expand_only_without_membrane))
    def nuclear_mask_flag = segment_nuclei ? "--nuclear-mask ${nuclear_mask}" : ''
    """
    expand_cells.py \\
        --whole-cell-mask ${whole_cell_mask} \\
        ${nuclear_mask_flag} \\
        --output ${prefix}_whole-cell_expanded.tiff \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_whole-cell_expanded.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
    END_VERSIONS
    """
}
