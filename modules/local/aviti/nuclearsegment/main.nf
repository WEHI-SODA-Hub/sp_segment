process AVITINUCLEARSEGMENT {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/python_numpy_pandas_scikit-image_pruned:4cec16761766b006'
        : 'community.wave.seqera.io/library/python_numpy_pandas_scikit-image_pruned:4cec16761766b006'}"

    input:
    tuple val(meta), path(nucleus_tif)
    path model_path

    output:
    // Named after the source tile, not meta.id: this is exactly the AVITI
    // viewer filename (`<tile>_Nuclear.tif`) expected under `Well<well>/`.
    tuple val(meta), path("${meta.tile}_Nuclear.tif")      , emit: nuclear_mask
    // Instance-labeled (uint16) counterpart, pre-binarize -- not part of the
    // Elembio viewer/cells2stats contract, purely an internal artifact for
    // AVITISTITCHWELL + CELLMEASUREMENT, which need per-nucleus instance IDs
    // rather than a 0/1 presence mask.
    tuple val(meta), path("${meta.tile}_Nuclear_label.tif"), emit: nuclear_label_mask
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Diameter comes solely from ext.args (aviti_nuclear_diameter): the
    // nuclear model/diameter is one fixed choice for the whole run, not a
    // per-row property, so meta.cell_diameter (a per-row samplesheet
    // override intended for the whole-cell/membrane path) does not apply here.
    """
    aviti_cellpose3_segment.py --mode nuclear \\
        --nucleus-tif ${nucleus_tif} \\
        --output ${meta.tile}_Nuclear.tif \\
        --output-label ${meta.tile}_Nuclear_label.tif \\
        --model-path ${model_path} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.tile}_Nuclear.tif
    touch ${meta.tile}_Nuclear_label.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(python3 -c "from importlib.metadata import version; print(version('cellpose'))")
        torch: \$(python3 -c "import torch; print(torch.__version__)")
    END_VERSIONS
    """
}
