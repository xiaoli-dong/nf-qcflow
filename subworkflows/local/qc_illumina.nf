include {
    FASTQC as FASTQC_INPUT ;
    FASTQC as FASTQC_TRIMMED_WITH_BBDUK ;
    FASTQC as FASTQC_TRIMMED_WITH_FASTP ;
    FASTQC as FASTQC_DEHOST
} from '../../modules/nf-core/fastqc/main'

include { BBMAP_BBDUK } from '../../modules/nf-core/bbmap/bbduk/main'

include { FASTP } from '../../modules/nf-core/fastp/main.nf'

include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT ;
    SEQKIT_STATS as SEQKIT_STATS_TRIMMED_WITH_BBDUK ;
    SEQKIT_STATS as SEQKIT_STATS_TRIMMED_WITH_FASTP ;
    SEQKIT_STATS as SEQKIT_STATS_DEHOST
} from '../../modules/local/seqkit/stats/main.nf'

include { HTML_COPYDIR } from '../../modules/local/html/copydir/main.nf'

include {
    HTML_DATATABLE_CSV2JSON as HTML_DATATABLE_CSV2JSON_QCREPORT ;
    HTML_DATATABLE_CSV2JSON as HTML_DATATABLE_CSV2JSON_TOPMATCHES
} from '../../modules/local/html/datatable/csv2json/main.nf'


include {
    CSVTK_CONCAT as CSVTK_CONCAT_REPORT ;
    CSVTK_CONCAT as CSVTK_CONCAT_TOPMATCHES
} from '../../modules/nf-core/csvtk/concat/main.nf'

include { HOSTILE_CLEAN as HOSTILE_CLEAN_ILLUMINA } from '../../modules/local/hostile/clean/main.nf'

include { KRAKEN2_KRAKEN2 } from '../../modules/local/kraken2/kraken2/main.nf'
include { BRACKEN_BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { BRACKEN_GETTOPMATCHES } from '../../modules/local/bracken/gettopmatches/main'

include { REPORT_QCSUMMARY } from '../../modules/local/report/qcsummary/main'
workflow QC_ILLUMINA {
    take:
    short_reads

    main:

    ch_versions = Channel.empty()

    trimmed_reads = short_reads
    qc_reads = short_reads
    ch_input_to_qc_summary = Channel.empty()

    short_reads = short_reads.filter { meta, reads ->
        def readsList = reads instanceof List ? reads : [reads]
        readsList.size() > 0 && readsList.every { it != null && it.exists() && it.size() > 0 }
    }
    
    //short_reads.filter { it[1][0].size() > 0 && it[1][0].countFastq() > 0 }
    FASTQC_INPUT(short_reads)
    ch_versions = ch_versions.mix(FASTQC_INPUT.out.versions.first())


    SEQKIT_STATS_INPUT(short_reads)
    ch_input_to_qc_summary = SEQKIT_STATS_INPUT.out.stats
    ch_versions = ch_versions.mix(SEQKIT_STATS_INPUT.out.versions.first())


    qc_stats = Channel.empty()

    //default
    if (params.illumina_reads_qc_tool == 'bbduk') {
        BBMAP_BBDUK(short_reads, [])
        ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        //get rid of zero size contig file and avoid the downstream crash
        BBMAP_BBDUK.out.reads
            .filter { it[1][0].size() > 0 && it[1][0].countFastq() > 0 }
            .set { trimmed_reads }

        FASTQC_TRIMMED_WITH_BBDUK(trimmed_reads)
        SEQKIT_STATS_TRIMMED_WITH_BBDUK(trimmed_reads)


        qc_reads = trimmed_reads
        qc_stats = SEQKIT_STATS_TRIMMED_WITH_BBDUK.out.stats
        ch_input_to_qc_summary = ch_input_to_qc_summary.join(qc_stats)
    }
    else if (params.illumina_reads_qc_tool == 'fastp') {
        save_trimmed_fail = false
        discard_trimmed_pass = false
        save_merged = false
        FASTP(short_reads, [], discard_trimmed_pass, save_trimmed_fail, save_merged)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        //get rid of zero size contig file and avoid the downstream crash

        FASTP.out.reads
            .filter { meta, reads ->
                def readsList = reads instanceof List ? reads : [reads]
                readsList.size() > 0 && readsList.every { it != null && it.exists() && it.size() > 0 }
            }
            .set { trimmed_reads }

        //qc_reads = FASTP.out.reads
        FASTQC_TRIMMED_WITH_FASTP(trimmed_reads)
        SEQKIT_STATS_TRIMMED_WITH_FASTP(trimmed_reads)


        qc_reads = trimmed_reads
        qc_stats = SEQKIT_STATS_TRIMMED_WITH_FASTP.out.stats
        ch_input_to_qc_summary = ch_input_to_qc_summary.join(qc_stats)
    }


    //FASTQC_QC(qc_reads)
    if (!params.skip_dehost) {

        HOSTILE_CLEAN_ILLUMINA(trimmed_reads, [params.hostile_ref_name_illumina, params.hostile_ref_dir])
        ch_versions = ch_versions.mix(HOSTILE_CLEAN_ILLUMINA.out.versions.first())
        HOSTILE_CLEAN_ILLUMINA.out.fastq
            .filter { meta, reads ->
                def readsList = reads instanceof List ? reads : [reads]
                readsList.size() > 0 && readsList.every { it != null && it.exists() && it.size() > 0 }
            }
            .set { qc_reads }

        FASTQC_DEHOST(qc_reads)
        SEQKIT_STATS_DEHOST(qc_reads)

        //qc_reads = HOSTILE_ILLUMINA.out.reads
        qc_stats = SEQKIT_STATS_DEHOST.out.stats
        ch_input_to_qc_summary = ch_input_to_qc_summary.join(qc_stats)
    }
    //ch_input_to_qc_summary.view()

    KRAKEN2_KRAKEN2(
        qc_reads,
        Channel.value(file(params.kraken2_db)),
        false,
        true,
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, params.kraken2_db)
    BRACKEN_GETTOPMATCHES(BRACKEN_BRACKEN.out.reports)
    CSVTK_CONCAT_TOPMATCHES(
        BRACKEN_GETTOPMATCHES.out.csv.map { meta, csv -> csv }.collect().map { csvs ->
            tuple([id: "reads.topmatches"], csvs)
        },
        'csv',
        'csv',
    )

    ch_to_combine_bracken_report = BRACKEN_BRACKEN.out.reports
        .map { meta, report ->
            report
        }
        .collect()
        .map { reports ->
            tuple([id: "reads_illumina_bracken_report"], reports)
        }
    BRACKEN_COMBINEBRACKENOUTPUTS(ch_to_combine_bracken_report)
    ch_versions = ch_versions.mix(BRACKEN_COMBINEBRACKENOUTPUTS.out.versions)

    REPORT_QCSUMMARY(ch_input_to_qc_summary)
    ch_versions = ch_versions.mix(REPORT_QCSUMMARY.out.versions)

    CSVTK_CONCAT_REPORT(
        REPORT_QCSUMMARY.out.csv.map { it ->
            it[1]
        }.collect().map { files ->
            tuple([id: "reads.qc_report"], files)
        },
        'csv',
        'csv',
    )
    ch_versions = ch_versions.mix(CSVTK_CONCAT_REPORT.out.versions)

    ch_html_report_template = Channel.fromPath("${projectDir}/assets/html_report_template", checkIfExists: true)

    HTML_COPYDIR(ch_html_report_template, "illumina")

    HTML_DATATABLE_CSV2JSON_QCREPORT(CSVTK_CONCAT_REPORT.out.csv)
    ch_versions = ch_versions.mix(HTML_DATATABLE_CSV2JSON_QCREPORT.out.versions)

    HTML_DATATABLE_CSV2JSON_TOPMATCHES(CSVTK_CONCAT_TOPMATCHES.out.csv)
    ch_versions = ch_versions.mix(HTML_DATATABLE_CSV2JSON_TOPMATCHES.out.versions)

    emit:
    input_stats = SEQKIT_STATS_INPUT.out.stats
    qc_reads
    qc_stats
    versions = ch_versions
}
