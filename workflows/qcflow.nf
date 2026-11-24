/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
include { paramsSummaryMap } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_qcflow_pipeline'
include { QC_SHORT } from '../subworkflows/local/qc_short'
include { QC_LONG } from '../subworkflows/local/qc_long'


include {
    COPY_HTML_DIR ;
    COPY_VERSIONS_TO_HTMLDIR ;
    CSV2JSON_FOR_HTML_DATATABLE as CSV2JSON_FOR_HTML_DATATABLE_QCREPORT ;
    CSV2JSON_FOR_HTML_DATATABLE as CSV2JSON_FOR_HTML_DATATABLE_TOPMATCHES
} from '../modules/local/html.nf'

include {
    CSVTK_CONCAT as CSVTK_CONCAT_COMBINED_QCREPORTS ;
    CSVTK_CONCAT as CSVTK_CONCAT_COMBINED_TOPMATCHES
} from '../modules/nf-core/csvtk/concat/main'

include { PUBLISH_SHORTREADS } from '../modules/local/publish/shortreads'
include { PUBLISH_LONGREADS } from '../modules/local/publish/longreads'
include { PUBLISH_SAMPLESHEET } from '../modules/local/publish/samplesheet/main.nf'
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow QCFLOW {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    def ch_versions = channel.empty()
    //ch_samplesheet.view { "SAMPLESHEET: ${it}" }

    /*
     * Split by sequencing platform
     */

    ch_illumina = ch_samplesheet.filter { meta, _fastqs ->
        meta.platform == 'illumina'
    }

    ch_nanopore = ch_samplesheet.filter { meta, _fastqs ->
        meta.platform == 'nanopore'
    }

    /*
     * Platform-specific trimming
     */

    QC_SHORT(ch_illumina)
    QC_LONG(ch_nanopore)

    short_reads_collected = QC_SHORT.out.qcreads
        .collect()
        .ifEmpty([])
        .map { it.collate(2) }
        //.first()

     //short_reads_collected = PUBLISH_SHORTREADS.out.reads.collect().ifEmpty([]).collate(2)
    long_reads_collected = QC_LONG.out.qcreads
        .collect()
        .ifEmpty([])
        .map { it.collate(2) }

    PUBLISH_SAMPLESHEET(
        short_reads_collected,
        long_reads_collected
    )

    ch_html_report_template = channel.fromPath("${projectDir}/assets/html_report_template", checkIfExists: true)
    COPY_HTML_DIR(ch_html_report_template)

    // 1. Combine the two channels into one
    ch_combined_reports = QC_SHORT.out.report_qc.mix(QC_LONG.out.report_qc).collect()
    //ch_combined_reports.view { "COMBINED QC REPORTS: ${it}" }

    // 3. Pass the list to CSVTK_CONCAT
    CSVTK_CONCAT_COMBINED_QCREPORTS(
        ch_combined_reports.map { files -> tuple([id: "reads.qcreport"], files) },
        "csv",
        "csv",
    )

    CSV2JSON_FOR_HTML_DATATABLE_QCREPORT(CSVTK_CONCAT_COMBINED_QCREPORTS.out.csv)

    // 1. Combine the two channels into one
    ch_combined_topmatches = QC_SHORT.out.report_topmatches
        .map { meta, csv -> csv }
        .mix(
            QC_LONG.out.report_topmatches.map { meta, csv -> csv }
        )
        .collect()
    //ch_combined_topmatches.view { "COMBINED TOPMATCHES REPORTS: ${it}" }

    // 3. Pass the list to CSVTK_CONCAT
    CSVTK_CONCAT_COMBINED_TOPMATCHES(
        ch_combined_topmatches.map { files -> tuple([id: "reads.topmatches"], files) },
        "csv",
        "csv",
    )
    CSV2JSON_FOR_HTML_DATATABLE_TOPMATCHES(CSVTK_CONCAT_COMBINED_TOPMATCHES.out.csv)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'qcflow_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    COPY_VERSIONS_TO_HTMLDIR(
        ch_collated_versions
    )

    emit:
    versions = ch_collated_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
