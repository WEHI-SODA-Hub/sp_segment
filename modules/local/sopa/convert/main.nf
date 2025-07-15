process SOPA_CONVERT {
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
