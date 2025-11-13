include {
    NANOPLOT as NANOPLOT_INPUT ;
    NANOPLOT as NANOPLOT_QC
} from '../../modules/local/nanoplot/main'

include { PORECHOP_PORECHOP } from '../../modules/local/porechop/porechop/main.nf'
include { CHOPPER } from '../../modules/local/chopper/main.nf'
include { FASTPLONG } from '../../modules/local/fastplong/main.nf'
include { HOSTILE_CLEAN as HOSTILE_CLEAN_NANOPORE } from '../../modules/local/hostile/clean/main'
include { DEACON_FILTER as DEACON_FILTER_NANOPORE } from '../../modules/local/deacon/filter/main'

include {
    SEQKIT_STATS as SEQKIT_STATS_INPUT ;
    SEQKIT_STATS as SEQKIT_STATS_PORECHOP ;
    SEQKIT_STATS as SEQKIT_STATS_CHOPPER ;
    SEQKIT_STATS as SEQKIT_STATS_FASTPLONG ;
    SEQKIT_STATS as SEQKIT_STATS_DEHOST
} from '../../modules/local/seqkit/stats/main'

include { KRAKEN2_KRAKEN2 } from '../../modules/local/kraken2/kraken2/main.nf'
include { BRACKEN_BRACKEN } from '../../modules/local/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../../modules/local/bracken/combinebrackenoutputs/main'
include { BRACKEN_GETTOPMATCHES } from '../../modules/local/bracken/gettopmatches/main'


include { REPORT_QCSUMMARY } from '../../modules/local/report/qcsummary/main'
include {
    CSVTK_CONCAT as CSVTK_CONCAT_REPORT ;
    CSVTK_CONCAT as CSVTK_CONCAT_TOPMATCHES
} from '../../modules/local/csvtk/concat/main.nf'

//include { CSVTK_CONCAT } from '../../modules/local/csvtk/concat/main.nf'
include { HTML_COPYDIR } from '../../modules/local/html/copydir/main.nf'

include {
    HTML_DATATABLE_CSV2JSON as HTML_DATATABLE_CSV2JSON_QCREPORT ;
    HTML_DATATABLE_CSV2JSON as HTML_DATATABLE_CSV2JSON_TOPMATCHES
} from '../../modules/local/html/datatable/csv2json/main.nf'



workflow QC_NANOPORE {
    take:
    long_reads

    main:

    ch_versions = channel.empty()
    ch_all_stats = channel.empty()
    qc_reads =  long_reads

    reads = long_reads.filter { meta, reads ->
        def readsList = reads instanceof List ? reads : [reads]
        readsList.size() > 0 && readsList.every { it != null && it.exists() && it.size() > 0 }
    }

    NANOPLOT_INPUT(reads)
    ch_versions = ch_versions.mix(NANOPLOT_INPUT.out.versions.first())

    SEQKIT_STATS_INPUT(reads)
    ch_all_stats = SEQKIT_STATS_INPUT.out.stats
    ch_versions = ch_versions.mix(SEQKIT_STATS_INPUT.out.versions.first())

    // QC
    if (params.long_qc_tool == 'fastplong') {
        discard_trimmed_pass = false
        save_trimmed_fail = false

        if (params.contaminants_fa) {
            FASTPLONG(reads, params.contaminants_fa, discard_trimmed_pass, save_trimmed_fail)

            ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash

            FASTPLONG.out.reads
                .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
                .set { trimmed_reads }
        }
        else {
            FASTPLONG(reads, [], discard_trimmed_pass, save_trimmed_fail)
            ch_versions = ch_versions.mix(FASTPLONG.out.versions.first())
            //get rid of zero size contig file and avoid the downstream crash

            FASTPLONG.out.reads
                .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
                .set { trimmed_reads }
        }
        NANOPLOT_QC(trimmed_reads)
        ch_versions = ch_versions.mix(NANOPLOT_QC.out.versions.first())

        SEQKIT_STATS_FASTPLONG(trimmed_reads)
        qc_reads = trimmed_reads
        qc_stats = SEQKIT_STATS_FASTPLONG.out.stats
        ch_all_stats = ch_all_stats.join(qc_stats)
    }
    else if(params.long_qc_tool == 'porechop+chopper') {
        PORECHOP_PORECHOP(reads)
        ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions.first())
        PORECHOP_PORECHOP.out.reads
            .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
            .set { nanopore_reads }
        SEQKIT_STATS_PORECHOP(nanopore_reads)
        ch_versions = ch_versions.mix(SEQKIT_STATS_PORECHOP.out.versions.first())

        qc_stats = SEQKIT_STATS_PORECHOP.out.stats

        if (params.contaminants_fa) {
            CHOPPER(nanopore_reads, params.contaminants_fa)
            ch_versions = ch_versions.mix(CHOPPER.out.versions.first())

            CHOPPER.out.fastq
                .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
                .set { qc_reads }
        }
        else {
            CHOPPER(nanopore_reads, [])
            ch_versions = ch_versions.mix(CHOPPER.out.versions.first())

            CHOPPER.out.fastq
                .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
                .set { qc_reads }
        }
        NANOPLOT_QC(qc_reads)
        ch_versions = ch_versions.mix(NANOPLOT_QC.out.versions.first())

        SEQKIT_STATS_CHOPPER(qc_reads)
        qc_stats = SEQKIT_STATS_CHOPPER.out.stats
        ch_all_stats = ch_all_stats.join(qc_stats)
    }

    if (!params.skip_dehost) {

        if (params.dehoster == 'hostile') {
            HOSTILE_CLEAN_NANOPORE(qc_reads, [params. hostile_refdb_long, params.hostile_ref_dir])
            ch_versions = ch_versions.mix(HOSTILE_CLEAN_NANOPORE.out.versions.first())

            HOSTILE_CLEAN_NANOPORE.out.fastq
                .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
                .set { qc_reads }

            SEQKIT_STATS_DEHOST(qc_reads)
            qc_stats = SEQKIT_STATS_DEHOST.out.stats
            ch_all_stats = ch_all_stats.join(qc_stats)
            ch_versions = ch_versions.mix(SEQKIT_STATS_DEHOST.out.versions.first())
        }
        else if (params.dehoster == 'deacon') {

            DEACON_FILTER_NANOPORE(
                qc_reads.map{meta, reads -> [meta, params.deacon_refdb, reads]}
            )
            ch_versions = ch_versions.mix(DEACON_FILTER_NANOPORE.out.versions.first())

           /*  DEACON_FILTER_NANOPORE.out.fastq_filtered
                .filter { meta, reads -> reads.size() > 0 && reads.countFastq() > 0 }
                .set { qc_reads }
 */
            qc_reads = DEACON_FILTER_NANOPORE.out.fastq_filtered.filter { meta, reads ->
                def readsList = reads instanceof List ? reads : [reads]
                readsList.size() > 0 && readsList.every { it != null && it.exists() && it.size() > 0 }
            }

            SEQKIT_STATS_DEHOST(qc_reads)
            qc_stats = SEQKIT_STATS_DEHOST.out.stats
            ch_all_stats = ch_all_stats.join(qc_stats)
            ch_versions = ch_versions.mix(SEQKIT_STATS_DEHOST.out.versions.first())
        }
    }

    //merged input, trimmed, and dehost stats into one columns
    // qc_summary.py
    //ch_all_stats.view()

    KRAKEN2_KRAKEN2(qc_reads, params.kraken2_db, false, true)
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
            tuple([id: "reads_nanopore_bracken_report"], reports)
        }
    BRACKEN_COMBINEBRACKENOUTPUTS(ch_to_combine_bracken_report)
    ch_versions = ch_versions.mix(BRACKEN_COMBINEBRACKENOUTPUTS.out.versions)

    REPORT_QCSUMMARY(ch_all_stats)

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

    ch_html_report_template = channel.fromPath("${projectDir}/assets/html_report_template", checkIfExists: true)

    HTML_COPYDIR(ch_html_report_template, "nanopore")

    HTML_DATATABLE_CSV2JSON_QCREPORT(CSVTK_CONCAT_REPORT.out.csv)
    ch_versions = ch_versions.mix(HTML_DATATABLE_CSV2JSON_QCREPORT.out.versions)

    HTML_DATATABLE_CSV2JSON_TOPMATCHES(CSVTK_CONCAT_TOPMATCHES.out.csv)
    ch_versions = ch_versions.mix(HTML_DATATABLE_CSV2JSON_TOPMATCHES.out.versions)

    emit:
    qc_reads
    versions = ch_versions
}
