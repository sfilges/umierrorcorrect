# UMIErrorCorrect

[![PyPI version](https://badge.fury.io/py/umierrorcorrect.svg)](https://badge.fury.io/py/umierrorcorrect)
[![Python Versions](https://img.shields.io/pypi/pyversions/umierrorcorrect.svg)](https://pypi.org/project/umierrorcorrect/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Pipeline for analyzing barcoded amplicon sequencing data with Unique Molecular Identifiers (UMI). This is essentially a complete overhaul of the original [UMIErrorCorrect](https://github.com/stahlberggroup/umierrorcorrect) published in [Clinical Chemistry](https://doi.org/10.1093/clinchem/hvac136).

The following major changes have been made:

- Added fastp as the default preprocessor including read trimming, quality filtering, merging and error correction of paired-end reads etc.
- Included quality control with fastqc and report merging with multiqc
- Added batch processing of fastq files and samplesheet support
- A complete refactoring of the code base, including removal of legacy and ununsed code, modernization of the codebase, and improved documentation
- A new command line interface using typer + rich
- Added a comprehensive testing suite of >100 tests

## Installation

### Using pip

```bash
pip install umierrorcorrect
```

### Using uv (recommended)

```bash
uv pip install umierrorcorrect
```

### Using Docker

```bash
docker pull ghcr.io/sfilges/umierrorcorrect:latest
# Or build locally
docker build -t umierrorcorrect docker/
```

See the [Docker documentation](doc/docker.md) for more details.

### Verify installation

```bash
umierrorcorrect --help
```

## Tutorial

This tutorial shows an example on how to use `umierrorcorrect`.

### Step 1 - Download the example files

1. [Download the example files](https://guse.sharefile.eu/d-s878525e0e6234cd69de7b4e4cfa4c113)
2. Create a folder named `tutorial` and move the files there.

### Step 2 - Reference genome

The `umierrorcorrect` pipeline runs `bwa mem` for aligning the reads to the reference genome. In this tutorial the data contains reads from chr3 of the Human genome.

If you don't have a version of the Human reference genome, you can download the sequence in fasta format for chromosome 3 (or for the full data set) from the [UCSC genome browser](http://hgdownload.cse.ucsc.edu/downloads.html#human).

After downloading the reference fasta file, unzip the file and index it using `bwa`:

```bash
gunzip chr3.fa.gz
bwa index chr3.fa
```

### Step 3 - Run UMI error correct

The data in the file `test_1.fastq.gz` is single end data generated with the SimSenSeq protocol described in [Ståhlberg et al. (2017)](https://doi.org/10.1038/nprot.2017.006).

According to the protocol, the UMI length (`-ul`) should be 12 and spacer length (`-sl`) is 16.

Run the pipeline with the following command:

```bash
umierrorcorrect run -o tutorial_output -r1 test_1.fastq.gz \
    -r chr3.fa -bed target_regions_hg38.bed -ul 12 -sl 16 -t 1
```

**Parameters matches:**

- `-o`: Output directory (created if it doesn't exist).
- `-r1`: Path to the fastq file with sequencing data.
- `-r`: Path to the reference genome fasta file.
- `-bed`: Path to the bedfile for annotation (target regions).
- `-ul`: UMI sequence length (12).
- `-sl`: Spacer length (16).
- `-t`: Number of CPU threads.

### Step 4 - Filter the cons file

The pipeline should have created output files in the folder `tutorial_output`, including consensus files and variant calls.

## Dependencies

UMIErrorCorrect requires Python 3.9+ and the following:

**Python libraries** (installed automatically):

- pysam (≥0.8.4)
- scipy
- matplotlib
- loguru
- typer

**External programs** (must be in PATH):

- `bwa` - for read mapping
- `gzip` or `pigz` - for compression

### Reference genome setup

The pipeline uses `bwa` for mapping, so you need an indexed reference genome:

```bash
bwa index -a bwtsw reference.fa
```

## Usage

### Full pipeline

```bash
umierrorcorrect run -r1 read1.fastq.gz -r2 read2.fastq.gz \
    -ul <umi_length> -sl <spacer_length> \
    -r reference.fa -o output_directory
```

### Pipeline steps

The `umierrorcorrect run` command performs:

1. **Preprocessing** - Extract UMI from reads and add to header
2. **Mapping** - Align reads to reference genome with BWA
3. **UMI clustering** - Group reads by UMI similarity
4. **Error correction** - Generate consensus for each UMI cluster
5. **Statistics** - Create consensus output file with per-position counts
6. **Variant calling** - Call variants from consensus data

### Running individual steps

Each step can be run independently using subcommands:

```bash
umierrorcorrect preprocess --help
umierrorcorrect mapping --help
umierrorcorrect consensus --help
umierrorcorrect stats --help
umierrorcorrect variants --help
umierrorcorrect filter-bam --help
umierrorcorrect filter-cons --help
```

## Documentation

- [Tutorial](https://github.com/stahlberggroup/umierrorcorrect/wiki/Tutorial)
- [UMI definition options](https://github.com/stahlberggroup/umierrorcorrect/wiki/UMI-definition-options)

## Reference

> Osterlund T., Filges S., Johansson G., Stahlberg A. *UMIErrorCorrect and UMIAnalyzer: Software for Consensus Read Generation, Error Correction, and Visualization Using Unique Molecular Identifiers*, Clinical Chemistry, 2022, hvac136
