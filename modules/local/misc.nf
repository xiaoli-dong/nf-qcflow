process SEQKIT_STATS_SUMMARY {
    tag "$meta.id"

    label 'process_low'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python%3A3.12'
        : 'quay.io/biocontainers/python:3.12'}"



    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqstats_summary.py \\
        --input $stats \\
        --output ${prefix}_summary.tsv
    """
}
process READ_QCREPORT {
    tag "read_qc_report.${platform}"
    label 'process_low'

    conda "conda-forge::pandas=2.2.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:2.2.1'
        : 'quay.io/biocontainers/pandas:2.2.1'}"

    input:
    val platform
    path input_summary
    path trimmed_summary
    path dehosted_summary

    output:
    path "${prefix}_${platform}.csv", emit: csv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions

    script:
    prefix = task.ext.prefix ?: "read.qcreport"
    """
    read_qc_report.py \\
        --input ${input_summary} \\
        --trimmed ${trimmed_summary} \\
        --dehosted ${dehosted_summary} \\
        --output ${prefix}_${platform}.csv
    """
}

