/*
 * This module uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
process SOPA_RESOLVECELLPOSE {
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.0.7'
        : 'docker.io/quentinblampey/sopa:2.0.7'}"

    input:
    tuple val(meta), path(zarr), val(cellpose_parquet)

    output:
    tuple val(meta), path("*.zarr/shapes/cellpose_boundaries/*.parquet"), emit: cellpose_boundaries

    script:
    """
    sopa resolve cellpose ${zarr}
    """
}
