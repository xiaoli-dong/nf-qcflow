# nf-qcflow
nf-qcflow is a bioinformatics pipeline that can be used to do NGS sequence quality control, host removal, initial contamination checking, and generate sequence tabular format statistical reports and html reports. The pipeline works with both Illumina and Nanopore data

## Pipeline summary

nf-qcflow supports both short and long reads:


- Short reads quality check and quality control
  - Short Illumina reads quality checks ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
  - Short read quality control ([BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) | [fastp](https://github.com/OpenGene/fastp))
  - dehost ([hostile](https://github.com/bede/hostile))
  - Taxonomic assignment and contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
  - Short read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
  - Tabular and html reports
- Long reads quality check and quality control
  - Nanopore long read quality checks ([NanoPlot](https://github.com/wdecoster/NanoPlot))
  - Nanopore long read adapter trimming, quality and length filter (porechop+chopper | fastplong)
     - Porechop + Chopper 
        - Nanopore long reads adapter removal ([Porechop](https://github.com/rrwick/Porechop))
        - Nanopore long read quality and length filter ([chopper](https://github.com/wdecoster/chopper))
    - Ultrafast preprocessing and quality control for long reads ([fastplong](https://github.com/OpenGene/fastplong))
  - dehost ([hostile](https://github.com/bede/hostile))
  - Taxonomic assignment and contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
  - Nanopore long read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
  - Tabular and html reports

## Quick start

>If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with -profile test before running the workflow on actual data.

### Check workflow options
You can clone or download the nf-qcflow from github to local computer or you can directly run the pipeline from github. To check the pipeline command line options:

```{r df-drop-ok, class.source="bg-success"}
# running directly from github without downloading or cloning
nextflow run xiaoli-dong/nf-qcflow -r revision_number(e.g:04b8745) --help
```
### Prepare required samplesheet input
The nf-qcflow pipeline requires user to provide a csv format samplesheet, which contains the sequenence information for each sample, as input. See below for what the samplesheet looks like:

`samplesheet.csv` for paired-end data:

```csv
sample,fastq_1,fastq_2
sample1,./fastq/3_S3_L001_R1_001.fastq.gz,./fastq/3_S3_L001_R2_001.fastq.gz
sample2,./fastq/78_S30_L001_R1.fastq.gz,./fastq/78_S30_L001_R2.fastq.gz
measle,./fastq/47_S47_L001_R1_001.fastq.gz,./fastq/47_S47_L001_R2_001.fastq.gz
flua,./fastq/30_S30_L001_R1_001.fastq.gz,./fastq/30_S30_L001_R2_001.fastq.gz
flub,./fastq/315_S122_L001_R1_001.fastq.gz,./fastq/315_S122_L001_R2_001.fastq.gz
```

`samplesheet.csv` for single-end data:

```csv
sample,fastq_1,fastq_2
sample1,./fastq/barcode01.fastq.gz,
sample2,./fastq/barcode02.fastq.gz,
```

The csv format samplesheet has three columns:
* The first row of the csv file is the header describing the columns
* Each row represents a unique sample to be processed, the first colum is the unique sample id
* Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

### Run the pipeline:
Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash

# Example command to run the pipeline from local download for illumina data
nextflow run nf-qcflow/main.nf \
  -profile singularity \
  --input samplesheet.csv \
  --platform illumina \
  --outdir results \
  -resume

# Example command to run the pipeline from local download for nanopore data
nextflow run nf-qcflow/main.nf \
  -profile singularity \
  --input samplesheet.csv \
  --platform nanopore \
  --outdir results \
  -resume

```

>* Notes: Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits
xiaoli-dong/qcflow was originally written by Xiaoli Dong.
