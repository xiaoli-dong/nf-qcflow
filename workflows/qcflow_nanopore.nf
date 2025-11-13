/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_qcflow_pipeline'

include { QC_NANOPORE } from '../subworkflows/local/qc_nanopore'

include { PUBLISH_LONGREADS } from '../modules/local/publish/longreads'
include { PUBLISH_SAMPLESHEET } from '../modules/local/publish/samplesheet/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QCFLOW_NANOPORE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_samplesheet.view()

    QC_NANOPORE(
            ch_samplesheet,
    )
    //QC_NANOPORE.out.versions.view()
    ch_versions = ch_versions.mix(QC_NANOPORE.out.versions)
    // Publish reads (only if channels have data)
    PUBLISH_LONGREADS(QC_NANOPORE.out.qc_reads)
    long_reads_collected = PUBLISH_LONGREADS.out.reads
        .collect()
        .ifEmpty([])
        .map { it.collate(2) }
        //.first()

    // Generate samplesheet
    PUBLISH_SAMPLESHEET(
        Channel.value([]),
        long_reads_collected
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf-qcflow_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )



    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
