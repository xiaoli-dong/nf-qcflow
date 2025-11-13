// Process to generate sample sheet
process PUBLISH_SAMPLESHEET {
    publishDir "${params.outdir}/samplesheets", mode: 'copy'

    input:
    val short_reads_list
    val long_reads_list

    output:
    path "*samplesheet.csv"

    exec:
    def short_reads_map = [:]
    def long_reads_map = [:]
    def prefix = task.ext.prefix ?: ""
    outdir_abs_path = file(params.outdir).toAbsolutePath().toString() + "/qcreads"
    // Parse short reads
    if (short_reads_list && short_reads_list instanceof List) {
        short_reads_list.each { item ->
            if (item instanceof List && item.size() >= 2) {
                def meta = item[0]
                def files = item[1]

                if (meta?.id && files) {
                    def fastq_1 = files instanceof List ?
                        new File(files[0].toString()).getName() :
                        new File(files.toString()).getName()
                    def fastq_2 = (files instanceof List && files.size() > 1) ?
                        new File(files[1].toString()).getName() : "NA"

                    short_reads_map[meta.id] = [
                        fastq_1: fastq_1,
                        fastq_2: fastq_2
                    ]
                }
            }
        }
    }

    // Parse long reads
    if (long_reads_list && long_reads_list instanceof List) {
        long_reads_list.each { item ->
            if (item instanceof List && item.size() >= 2) {
                def meta = item[0]
                def file = item[1]

                if (meta?.id && file) {
                    def filename = file instanceof List ?
                        new File(file[0].toString()).getName() :
                        new File(file.toString()).getName()
                    def basecaller = meta.basecaller_mode ?: "NA"

                    long_reads_map[meta.id] = [
                        long_fastq: filename,
                        basecaller_mode: basecaller
                    ]
                }
            }
        }
    }

    // Combine samples with same meta.id
    def all_samples = (short_reads_map.keySet() + long_reads_map.keySet()).unique().sort()
    def rows = ["sample,fastq_1,fastq_2,long_fastq,basecaller_mode"]

    all_samples.each { sample_id ->
        def short_data = short_reads_map[sample_id]
        def long_data = long_reads_map[sample_id]

        def fastq_1 = short_data?.fastq_1 ?: "NA"
        def fastq_2 = short_data?.fastq_2 ?: "NA"
        def long_fastq = long_data?.long_fastq ?: "NA"
        def basecaller_mode = long_data?.basecaller_mode ?: "NA"

        def output_fastq_1 = (fastq_1 != 'NA') ? "${outdir_abs_path}/${fastq_1}" : 'NA'
        def output_fastq_2 = (fastq_2 != 'NA') ? "${outdir_abs_path}/${fastq_2}" : 'NA'
        def output_long_fastq = (long_fastq != 'NA') ? "${outdir_abs_path}/${long_fastq}" : 'NA'
        rows << "${sample_id},${output_fastq_1},${output_fastq_2},${output_long_fastq},${basecaller_mode}"
    }

    // Write samplesheet to output
    def samplesheet = task.workDir.resolve("${prefix}samplesheet.csv")
    samplesheet.text = rows.join('\n')
}
