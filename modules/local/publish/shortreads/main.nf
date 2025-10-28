// Publish QC reads to subdirectory
    // Process to publish short reads
process PUBLISH_SHORTREADS {
    tag "$meta.id"
    //publishDir "${params.outdir}/qc_reads/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), emit: reads

    script:
    """
    # Files are automatically published by publishDir
    """
}

