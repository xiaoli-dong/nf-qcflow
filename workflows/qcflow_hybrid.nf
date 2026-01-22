/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_qcflow_pipeline'

include { QC_NANOPORE } from '../subworkflows/local/qc_nanopore'
include { QC_ILLUMINA } from '../subworkflows/local/qc_illumina'

include { PUBLISH_SHORTREADS } from '../modules/local/publish/shortreads'
include { PUBLISH_LONGREADS } from '../modules/local/publish/longreads'
include { PUBLISH_SAMPLESHEET } from '../modules/local/publish/samplesheet/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QCFLOW_HYBRID {
    take:
    short_reads
    long_reads // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    //short_reads.view()
    //long_reads.view()


    QC_ILLUMINA(short_reads)
    ch_versions = ch_versions.mix(QC_ILLUMINA.out.versions)

    QC_NANOPORE(long_reads)
    ch_versions = ch_versions.mix(QC_NANOPORE.out.versions)

    PUBLISH_SHORTREADS(QC_ILLUMINA.out.qc_reads)
    PUBLISH_LONGREADS(QC_NANOPORE.out.qc_reads)

     //short_reads_collected = PUBLISH_SHORTREADS.out.reads.collect().ifEmpty([]).collate(2)
    short_reads_collected = PUBLISH_SHORTREADS.out.reads
        .collect()
        .ifEmpty([])
        .map { it.collate(2) }
        //.first()

     //short_reads_collected = PUBLISH_SHORTREADS.out.reads.collect().ifEmpty([]).collate(2)
    long_reads_collected = PUBLISH_LONGREADS.out.reads
        .collect()
        .ifEmpty([])
        .map { it.collate(2) }
        //.first()

    PUBLISH_SAMPLESHEET(
        short_reads_collected,
        long_reads_collected
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'software_versions.yml',
        sort: true,
        newLine: true,
    )

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
