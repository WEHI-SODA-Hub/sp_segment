process AVITIMEMBRANESEGMENT {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/python_numpy_pandas_scikit-image_pruned:4cec16761766b006'
        : 'community.wave.seqera.io/library/python_numpy_pandas_scikit-image_pruned:4cec16761766b006'}"

    input:
    // model_path is carried in the same tuple as meta/tiles, not passed as a
    // separate input: unlike the nuclear model (one shared file for the
    // whole run, broadcast via a value channel), this row's membrane model
    // varies per meta.membrane_model, so it is a genuinely per-tile value.
    tuple val(meta), path(nucleus_tif), path(membrane_tif), path(actin_tif), path(model_path)

    output:
    tuple val(meta), path("${meta.tile}_Cell.tif"), emit: cell_mask
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Row-level cell_diameter must come *after* args: both set --diameter,
    // and typer/click keeps the last occurrence of a repeated option, so
    // this order is what lets a per-row override actually win over the
    // global aviti_membrane_diameter in ext.args.
    def diameter_arg = (meta.cell_diameter != null && meta.cell_diameter.toFloat() > 0) ? "--diameter ${meta.cell_diameter}" : ''
    def actin_arg = (actin_tif.name != 'NO_FILE') ? "--actin-tif ${actin_tif}" : ''
    """
    aviti_cellpose3_segment.py --mode membrane \\
        --cell-tif ${membrane_tif} \\
        --nucleus-tif ${nucleus_tif} \\
        ${actin_arg} \\
        --output ${meta.tile}_Cell.tif \\
        --model-path ${model_path} \\
        ${args} \\
        ${diameter_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tile}_Cell.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """
}
