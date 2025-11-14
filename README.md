# nf-qcflow

A comprehensive Nextflow bioinformatics pipeline for next-generation sequencing (NGS) data quality processing. nf-qcflow performs quality control, host sequence removal, contamination screening, and generates detailed statistical reports in both tabular and HTML formats. The pipeline supports both Illumina short-read and Oxford Nanopore long-read sequencing platforms.

---

## Pipeline Summary

nf-qcflow supports both short and long reads:

### Short Reads Quality Check and Quality Control

- Short Illumina reads quality checks ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Short read quality control ([BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) | [fastp](https://github.com/OpenGene/fastp))
- Dehost ([hostile](https://github.com/bede/hostile) | [deacon](https://github.com/dcdanko/deacon))
- Taxonomic assignment and contamination check ([Kraken2](https://ccb.jhu.edu/software/kraken2/))
- Short read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
- Tabular and HTML reports

### Long Reads Quality Check and Quality Control

- Nanopore long read quality checks ([NanoPlot](https://github.com/wdecoster/NanoPlot))
- Nanopore long read adapter trimming, quality and length filtering:
  - **Porechop + Chopper**
    - Adapter removal ([Porechop](https://github.com/rrwick/Porechop))
    - Quality and length filter ([chopper](https://github.com/wdecoster/chopper))
  - **fastplong** - Ultrafast preprocessing and quality control ([fastplong](https://github.com/OpenGene/fastplong))
- Dehost ([hostile](https://github.com/bede/hostile) | [deacon](https://github.com/dcdanko/deacon))
- Taxonomic assignment and contamination check ([Kraken2](https://ccb.jhu.edu/software/kraken2/))
- Nanopore long read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
- Tabular and HTML reports

---

## Quick Start

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Check Workflow Options

You can clone or download nf-qcflow from GitHub to your local computer, or you can run the pipeline directly from GitHub. To check the pipeline command-line options:

```bash
# Running directly from GitHub without downloading or cloning
nextflow run xiaoli-dong/nf-qcflow -r <revision_number> --help

# Example with specific revision
nextflow run xiaoli-dong/nf-qcflow -r 04b8745 --help
```

### Prepare Required Samplesheet Input

The nf-qcflow pipeline requires a CSV format samplesheet containing sequence information for each sample. See below for what the samplesheet looks like:

**samplesheet.csv**

```csv
sample,fastq_1,fastq_2,long_fastq
sample_paired_short_long,shortreads_1.fastq.gz,shortreads_2.fastq.gz,longreads.fastq.gz
sample_single_short_long,shortreads.fastq,NA,longreads.fastq.gz
sample_only_long,NA,NA,longreads.fastq.gz
```

**Samplesheet Format Requirements:**

- The first row of the CSV file is the header describing the columns
- Each row represents a unique sample to be processed; the first column is the unique sample ID
- When information for a particular column is missing, fill the column with `NA`
- The `fastq_1` and `fastq_2` columns are reserved for supplying the short sequence files
- The `long_fastq` column is reserved for supplying the long sequence file

### Run the Pipeline

Now you can run the pipeline using:

**For Illumina data:**

```bash
# Example command to run the pipeline from local download
nextflow run nf-qcflow/main.nf \
  -profile singularity \
  --input samplesheet.csv \
  --platform illumina \
  --outdir results \
  -resume
```

**For Nanopore data:**

```bash
# Example command to run the pipeline from local download
nextflow run nf-qcflow/main.nf \
  -profile singularity \
  --input samplesheet.csv \
  --platform nanopore \
  --outdir results \
  -resume
```

> [!IMPORTANT]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration **except for parameters**. See [documentation](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files) for more details.

---

## nf-qcflow Command-Line Options

### Input/Output Options

| Parameter | Type | Description |
|-----------|------|-------------|
| `--input` | string | Path to a CSV file containing sample information |
| `--outdir` | string | Output directory for results (must be an absolute path on cloud storage) |

### General Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--platform` | string | `illumina` | Sequencing platform (`illumina` or `nanopore`) |
| `--skip_dehost` | boolean | `false` | Skip the dehosting step |
| `--dehoster` | string | `deacon` | Dehosting tool to use (`deacon` or `hostile`) |
| `--short_qc_tool` | string | `fastp` | Short read QC tool (`fastp` or `bbduk`) |
| `--long_qc_tool` | string | `fastplong` | Long read QC tool (`fastplong` or `porechop_chopper`) |
| `--publish_dir_mode` | string | `copy` | Method for publishing results |
| `--contaminants_fa` | string | - | Path to custom contaminants FASTA file |

### Database Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--kraken2_db` | string | `/your_path/kraken2/k2_standard_16gb_20250402` | Path to Kraken2 database |
| `--deacon_refdb` | string | `/your_path/deacon/v3/panhuman-1.k31w15.idx` | Path to deacon reference database |
| `--hostile_ref_dir` | string | `/your_path/hostile2` | Path to hostile reference directory |
| `--hostile_refdb_short` | string | `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` | Hostile reference for short reads |
| `--hostile_refdb_long` | string | `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401.mmi` | Hostile reference for long reads |

### Help Options

| Parameter | Description |
|-----------|-------------|
| `--help` | Show help for all top-level parameters (can specify a parameter for detailed help) |
| `--help_full` | Show help for all non-hidden parameters |
| `--show_hidden` | Show hidden parameters (use with `--help` or `--help_full`) |

---

## Credits

nf-qcflow was originally written by Xiaoli Dong.

## Support

For issues, questions, or feature requests, please [open an issue](https://github.com/xiaoli-dong/nf-qcflow/issues) on GitHub.
