# nf-qcflow

A comprehensive Nextflow bioinformatics pipeline for next-generation sequencing (NGS) data quality processing. nf-qcflow performs quality control, host sequence removal, contamination screening, and generates detailed statistical reports in both tabular and HTML formats. The pipeline supports both Illumina short-read and Oxford Nanopore long-read sequencing platforms.

---

## Pipeline Summary

nf-qcflow supports both short and long reads:

### Short Reads Quality Check and Quality Control

- Short Illumina reads quality checks ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Short read quality control ([BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) | [fastp](https://github.com/OpenGene/fastp))
- Dehost ([hostile](https://github.com/bede/hostile) | [deacon](https://github.com/dcdanko/deacon))
- Taxonomic assignment and contamination check ([Kraken2](https://ccb.jhu.edu/software/kraken2/) & [Bracken](https://ccb.jhu.edu/software/bracken/))
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
- Taxonomic assignment and contamination check ([Kraken2](https://ccb.jhu.edu/software/kraken2/) & [Bracken](https://ccb.jhu.edu/software/bracken/))
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
nextflow run xiaoli-dong/nf-qcflow -r main --help
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
nextflow run xiaoli-dong/nf-qcflow \
  -profile singularity \
  --input samplesheet.csv \
  --platform illumina \
  --outdir results \
  -resume
```

**For Nanopore data:**

```bash
nextflow run xiaoli-dong/nf-qcflow \
  -profile singularity \
  --input samplesheet.csv \
  --platform nanopore \
  --outdir results \
  -resume
```

> [!IMPORTANT]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration **except for parameters**. See [documentation](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files) for more details.

---

## Command-Line Options

### Input/Output Options

| Parameter | Type | Description |
|-----------|------|-------------|
| `--input` | string | Path to CSV file containing information about the samples |
| `--outdir` | string | Output directory where results will be saved (use absolute paths for cloud storage) |

### General Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--platform` | string | `illumina` | Sequencing platform (`illumina` or `nanopore`) |
| `--skip_dehost` | boolean | `false` | Skip the dehosting step |
| `--dehoster` | string | `deacon` | Dehosting tool to use (`deacon` or `hostile`) |
| `--short_qc_tool` | string | `fastp` | Short read QC tool (`fastp` or `bbduk`) |
| `--long_qc_tool` | string | `fastplong` | Long read QC tool (`fastplong` or `porechop_chopper`) |
| `--contaminants_fa` | string | - | Path to custom contaminants FASTA file |
| `--publish_dir_mode` | string | `copy` | Method for publishing results |

### Database Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--kraken2_db` | string | `/your_path_to/kraken2_db` | Path to Kraken2 database |
| `--deacon_refdb` | string | `/your_path_to/deacon_db` | Path to deacon reference database |
| `--hostile_ref_dir` | string | `/your_path_to/hostile2_db_dir` | Path to hostile reference directory |
| `--hostile_refdb_short` | string | `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401` | Hostile reference for short reads |
| `--hostile_refdb_long` | string | `human-t2t-hla.argos-bacteria-985_rs-viral-202401_ml-phage-202401.mmi` | Hostile reference for long reads |

### Help Options

| Parameter | Description |
|-----------|-------------|
| `--help` | Show help for all top-level parameters (can specify a parameter for detailed help) |
| `--help_full` | Show help for all non-hidden parameters |
| `--show_hidden` | Show hidden parameters (use with `--help` or `--help_full`) |

---

## Workflow Overview

```
Input Reads (FASTQ)
    ↓
Quality Assessment
  • Short reads: FastQC
  • Long reads: NanoPlot
    ↓
Quality Control & Trimming
  • Short reads: fastp or BBDuk
  • Long reads: fastplong or Porechop+Chopper
    ↓
Host Sequence Removal (deacon or hostile) [optional]
    ↓
Taxonomic Classification (Kraken2 & Bracken)
    ↓
Read Statistics (seqkit)
    ↓
Summary Reports (tabular & HTML)
    ↓
Final QC Outputs
```
<!--
---

## Output Structure

```
results/
├── fastqc/              # FastQC reports (Illumina)
├── nanoplot/            # NanoPlot reports (Nanopore)
├── qc/                  # Quality-controlled reads
├── dehost/              # Dehosted reads (if applicable)
├── kraken2/             # Kraken2 taxonomic classification results
├── bracken/             # Bracken abundance estimation results
├── stats/               # Read statistics
└── reports/             # Summary reports (HTML and tabular)
```
-->
---

## Database Requirements

nf-qcflow requires several databases to be downloaded and configured:

### Kraken2 Database
```bash
# Download a pre-built Kraken2 database
# Standard databases available at: https://benlangmead.github.io/aws-indexes/k2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20250402.tar.gz
tar -xvzf k2_standard_16gb_20250402.tar.gz

```

### Deacon Database
```bash
# Download deacon reference database
# Visit: https://github.com/dcdanko/deacon for instructions
```

### Hostile Database
```bash
# Hostile databases are downloaded automatically on first use
# Or download manually from: https://github.com/bede/hostile
```

Provide database paths using the respective parameters (e.g., `--kraken2_db`, `--deacon_refdb`, `--hostile_ref_dir`).

---

## Tool References

This pipeline uses the following tools:

- [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Quality control for high throughput sequence data
- [**fastp**](https://github.com/OpenGene/fastp) - Ultra-fast all-in-one FASTQ preprocessor
- [**BBDuk**](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) - Kmer-based quality and adapter trimming
- [**NanoPlot**](https://github.com/wdecoster/NanoPlot) - Plotting tools for long read sequencing data
- [**Porechop**](https://github.com/rrwick/Porechop) - Adapter trimmer for Oxford Nanopore reads
- [**Chopper**](https://github.com/wdecoster/chopper) - Quality and length filtering for long reads
- [**fastplong**](https://github.com/OpenGene/fastplong) - Ultra-fast preprocessing for long reads
- [**hostile**](https://github.com/bede/hostile) - Fast host removal tool
- [**deacon**](https://github.com/dcdanko/deacon) - Decontamination of sequencing data
- [**Kraken2**](https://ccb.jhu.edu/software/kraken2/) - Taxonomic sequence classification
- [**Bracken**](https://ccb.jhu.edu/software/bracken/) - Bayesian re-estimation of abundance with Kraken
- [**seqkit**](https://bioinf.shenwei.me/seqkit/) - Cross-platform toolkit for FASTA/Q file manipulation

---

## Citations

If you use nf-qcflow in your research, please cite the appropriate tools:

- **FastQC** - Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data.
- **fastp** - Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890.
- **NanoPlot** - De Coster, W., D'Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics*, 34(15), 2666-2669.
- **Chopper** - De Coster, W., & Rademakers, R. (2023). NanoPack2: population-scale evaluation of long-read sequencing data. *Bioinformatics*, 39(5), btad311.
- **hostile** - Bede, P. (2023). hostile: accurate decontamination of microbial sequences. *Journal of Open Source Software*.
- **Kraken2** - Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.
- **Bracken** - Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- **seqkit** - Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. *PLOS ONE*, 11(10), e0163962.

---

## Credits

nf-qcflow was originally written by Xiaoli Dong.

## Support

For issues, questions, or feature requests, please [open an issue](https://github.com/xiaoli-dong/nf-qcflow/issues) on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
