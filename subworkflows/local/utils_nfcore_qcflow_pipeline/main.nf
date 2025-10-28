//
// Subworkflow with functionality specific to the xiaoli-dong/qcflow pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { samplesheetToList } from 'plugin/nf-schema'
include { completionEmail } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to parse samplesheet
def parse_samplesheet(samplesheet_path) {
    def samples = []

    file(samplesheet_path)
        .readLines()
        .drop(1)
        .each { line ->
            def fields = line.split(',')
            def sample_id = fields[0]
            def fastq_1 = fields[1] == 'NA' ? null : fields[1]
            def fastq_2 = fields[2] == 'NA' ? null : fields[2]
            def long_fastq = fields[3] == 'NA' ? null : fields[3]
            def basecaller_mode = fields[4] == 'NA' ? null : fields[4]

            samples.add(
                [
                    sample: sample_id,
                    fastq_1: fastq_1,
                    fastq_2: fastq_2,
                    long_fastq: long_fastq,
                    basecaller_mode: basecaller_mode,
                ]
            )
        }

    return samples
}

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir //  string: The output directory where the results will be saved
    input //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    // Read and parse samplesheet with duplicate validation
    def seen_samples = [] as Set

    reads = Channel.fromPath(params.input)
        .splitText()
        .filter { line -> !line.trim().startsWith('#') && line.trim() != '' }
        .collect()
        .map { lines -> lines.join('\n') }
        .splitCsv(header: true, sep: ',')
        .map { row ->
            // Skip rows where sample column starts with #
            if (row.sample.startsWith('#')) {
                return null
            }

            // Check for duplicate sample IDs
            if (seen_samples.contains(row.sample)) {
                error("ERROR: Duplicate sample ID found in samplesheet: '${row.sample}'\n" +
                      "Each sample ID must be unique!\n" +
                      "Please check your input samplesheet: ${params.input}")
            }
            seen_samples.add(row.sample)

            // Return parsed row
            [
                sample: row.sample,
                fastq_1: row.fastq_1 == 'NA' ? null : row.fastq_1,
                fastq_2: row.fastq_2 == 'NA' ? null : row.fastq_2,
                long_fastq: row.long_fastq == 'NA' ? null : row.long_fastq,
                basecaller_mode: row.basecaller_mode == 'NA' ? null : row.basecaller_mode,
            ]
        }
        .filter { it != null }

    // Create channel for short reads only (where at least fastq_1 or fastq_2 is not null)
    short_reads = reads
        .filter { it.fastq_1 != null || it.fastq_2 != null }
        .map { row ->
            def meta = [
                id: row.sample,
                single_end: row.fastq_2 == null,
            ]
            def files = row.fastq_2 == null
                ? [file(row.fastq_1)]
                : [file(row.fastq_1), file(row.fastq_2)]
            tuple(meta, files)
        }

    // Create channel for long reads only (where long_fastq is not null)
    long_reads = reads
        .filter { it.long_fastq != null }
        .map { row ->
            def meta = [
                id: row.sample,
                single_end: true,
                basecaller_mode: row.basecaller_mode,
            ]
            tuple(meta, file(row.long_fastq))
        }

    // Convert samplesheet_short channel to CSV format
    short_reads
        .map { meta, reads ->
            // Handle both single-end and paired-end reads
            if (meta.single_end) {
                // Single-end: sample,fastq_1,NA
                return "${meta.id},${reads[0]},NA"
            }
            else {
                // Paired-end: sample,fastq_1,fastq_2
                return "${meta.id},${reads[0]},${reads[1]}"
            }
        }
        .collectFile(
            name: 'short_reads_samplesheet.csv',
            newLine: true,
            seed: 'sample,fastq_1,fastq_2',
            storeDir: "${params.outdir}/samplesheets",
        )
        .set { ch_short_reads_csv }

    // Log the generated samplesheet
    ch_short_reads_csv.view { "Generated samplesheet: ${it}" }

    emit:
    short_reads
    long_reads
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                [],
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastqs]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>"
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
