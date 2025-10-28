/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_qcflow_pipeline'
include { QC_ILLUMINA } from '../subworkflows/local/qc_illumina'
include { PUBLISH_SHORTREADS } from '../modules/local/publish/shortreads'
include { PUBLISH_SAMPLESHEET } from '../modules/local/publish/samplesheet'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QCFLOW_ILLUMINA {
    take:
    short_reads // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    //short_reads.view()
    QC_ILLUMINA(short_reads)
    ch_versions = ch_versions.mix(QC_ILLUMINA.out.versions)
    //QC_ILLUMINA.out.qc_reads.view()

    PUBLISH_SHORTREADS(QC_ILLUMINA.out.qc_reads)
    //short_reads_collected = PUBLISH_SHORTREADS.out.reads.collect().ifEmpty([]).collate(2)
    short_reads_collected = PUBLISH_SHORTREADS.out.reads
        .collect()
        .ifEmpty([])
        .map { it.collate(2) }
        //.first()

    short_reads_collected.view()
    // Generate samplesheet
    PUBLISH_SAMPLESHEET(
        short_reads_collected,
        Channel.value([]),
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions).collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'seqsweep_software_' + 'mqc_' + 'versions.yml',
        sort: true,
        newLine: true,
    )

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
