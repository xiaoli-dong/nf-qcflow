# nf-qcflow
nf-qcflow is a bioinformatics pipeline that can be used to do NGS sequence quality control, host removal, initial contamination checking, and generate sequence tabular format statistical reports and html reports. The pipeline works with both Illumina and Nanopore data

## Pipeline summary

nf-qcflow supports both short and long reads:

- Sequence quality check and quality control
  - Short reads
    - Short Illumina reads quality checks ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    - Short read quality control ([BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) | [fastp](https://github.com/OpenGene/fastp))
    - Short read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))
    - dehost ([hostile](https://github.com/bede/hostile))
    - Taxonomic assignment and contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
  - Long reads
    - Nanopore long read quality checks ([NanoPlot](https://github.com/wdecoster/NanoPlot))
    - Nanopore long read adapter trimming, quality and length filter (porechop+chopper | fastplong)
       - Porechop + Chopper 
          - Nanopore long reads adapter removal ([Porechop](https://github.com/rrwick/Porechop))
          - Nanopore long read quality and length filter ([chopper](https://github.com/wdecoster/chopper))
      - Ultrafast preprocessing and quality control for long reads ([fastplong](https://github.com/OpenGene/fastplong))
    - dehost ([hostile](https://github.com/bede/hostile))
    - Taxonomic assignment and contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
    - Nanopore long read statistics ([seqkit stats](https://bioinf.shenwei.me/seqkit/usage/#stats))

## Quick start

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
sample1,./fastq/3_S3_L001_R1_001.fastq.gz,./fastq/3_S3_L001_R2_001.fastq.gz
sample2,./fastq/78_S30_L001_R1.fastq.gz,./fastq/78_S30_L001_R2.fastq.gz
measle,./fastq/47_S47_L001_R1_001.fastq.gz,./fastq/47_S47_L001_R2_001.fastq.gz
flua,./fastq/30_S30_L001_R1_001.fastq.gz,./fastq/30_S30_L001_R2_001.fastq.gz
flub,./fastq/315_S122_L001_R1_001.fastq.gz,./fastq/315_S122_L001_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
#check command line options:
nextflow run nf-qcflow/main.nf --help

# example command to run the pipeline
nextflow run nf-qcflow/main.nf -profile singularity --input samplesheet.csv --platform illumina --outdir results
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

xiaoli-dong/qcflow was originally written by Xiaoli Dong.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use xiaoli-dong/qcflow for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
