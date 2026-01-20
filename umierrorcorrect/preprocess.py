#!/usr/bin/env python3
"""
UMI error correct, preprocess.py - remove UMI and append to the header of Fastq sequences.
================================

:Authors: Tobias Osterlund, Stefan Filges

Purpose
-------

Preprocess the fastq files by removing the unique molecular index and add it to the header of the fastq entry.

"""

import datetime
import subprocess
import tempfile
from pathlib import Path
from typing import Tuple

from umierrorcorrect.core.logging_config import get_logger, log_subprocess_stderr
from umierrorcorrect.core.read_fastq_records import read_fastq, read_fastq_paired_end
from umierrorcorrect.models.models import PreprocessConfig

logger = get_logger(__name__)


def generate_random_dir(tmpdir: str) -> str:
    """Generate a directory for storing temporary files, using a timestamp."""
    Path(tmpdir).mkdir(parents=True, exist_ok=True)
    prefix = f"r{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}_"
    return tempfile.mkdtemp(prefix=prefix, dir=str(tmpdir))


def run_unpigz(filename: str, tmpdir: str, num_threads: int, program: str) -> str:
    """Unzip fastq.gz files using parallel gzip (pigz) or gunzip."""
    input_path = Path(filename)
    outfilename = Path(tmpdir) / input_path.name.removesuffix(".gz")

    if program == "pigz":
        command = ["unpigz", "-p", str(num_threads), "-c", filename]
    else:
        command = ["gunzip", "-c", filename]

    with outfilename.open("w") as g:
        result = subprocess.run(command, stdout=g, stderr=subprocess.PIPE, check=True)
        log_subprocess_stderr(result.stderr, program)

    return str(outfilename)


def run_pigz(filename: str, num_threads: int, program: str) -> None:
    """Zip fastq files using parallel gzip (pigz) or gzip in-place."""
    if program == "pigz":
        command = ["pigz", "-p", str(num_threads), filename]
    else:
        command = ["gzip", filename]

    result = subprocess.run(command, capture_output=True, check=True)
    log_subprocess_stderr(result.stderr, program)


def run_cutadapt(
    r1file: str,
    r2file: str | None,
    output_path: Path,
    sample_name: str,
    adapter_sequence: str,
    mode: str,
) -> Tuple[str, str | None]:
    """Run cutadapt to trim adapters."""
    if adapter_sequence.lower() == "illumina":
        adapter = "AGATCGGAAGAGC"
    elif adapter_sequence.lower() == "nextera":
        adapter = "CTGTCTCTTATA"
    elif adapter_sequence.lower() == "small-rna":
        adapter = "ATGGAATTCTCG"
    else:
        adapter = adapter_sequence.upper()

    if mode == "single":
        outfilename = str(output_path / f"{sample_name}_trimmed.fastq")
        command = ["cutadapt", "-a", adapter, "-o", outfilename, "-O", "3", "-m", "20", r1file]
        outfile1 = outfilename
        outfile2 = None
    else:
        outfile1 = str(output_path / f"{sample_name}_R1_trimmed.fastq")
        outfile2 = str(output_path / f"{sample_name}_R2_trimmed.fastq")
        if r2file is None:
            raise ValueError("r2file cannot be None in paired mode")
        command = [
            "cutadapt",
            "-a",
            adapter,
            "-A",
            adapter,
            "-o",
            outfile1,
            "-p",
            outfile2,
            "-O",
            "3",
            "-m",
            "20",
            r1file,
            r2file,
        ]

    logger.info(f"Performing adapter trimming using cutadapt with adapter sequence {adapter}")
    result = subprocess.run(command, capture_output=True, check=True)
    log_subprocess_stderr(result.stderr, "cutadapt")

    # Clean up input files
    Path(r1file).unlink()
    if mode != "single" and r2file:
        Path(r2file).unlink()

    return outfile1, outfile2


def preprocess_se(infilename, outfilename, barcode_length, spacer_length):
    """Run the preprocessing for single end data (one fastq file)."""
    with Path(infilename).open() as f, Path(outfilename).open("w") as g:
        read_start = barcode_length + spacer_length
        nseqs = 0
        for name, seq, qual in read_fastq(f):
            nseqs += 1
            barcode = seq[:barcode_length]
            # g.write(name+':'+barcode+'\n'+rest+'\n'+qualname+'\n'+qual[12+11:]+'\n')
            parts = name.split()
            newname = ":".join([parts[0], barcode]) + " " + parts[-1]
            g.write("\n".join([newname, seq[read_start:], "+", qual[read_start:]]) + "\n")
    return nseqs


def preprocess_pe(r1file, r2file, outfile1, outfile2, barcode_length, spacer_length, dual_index):
    """Run the preprocessing for paired end data (two fastq files)."""
    read_start = barcode_length + spacer_length
    with (
        Path(r1file).open() as f1,
        Path(r2file).open() as f2,
        Path(outfile1).open("w") as g1,
        Path(outfile2).open("w") as g2,
    ):
        nseqs = 0
        for name1, seq1, qual1, name2, seq2, qual2 in read_fastq_paired_end(f1, f2):
            nseqs += 1
            if dual_index:
                barcode = seq1[:barcode_length] + seq2[:barcode_length]
            else:
                barcode = seq1[:barcode_length]
            parts1 = name1.split()
            parts2 = name2.split()
            newname1 = ":".join([parts1[0], barcode]) + " " + parts1[-1]
            newname2 = ":".join([parts2[0], barcode]) + " " + parts2[-1]
            g1.write("\n".join([newname1, seq1[read_start:], "+", qual1[read_start:]]) + "\n")
            if dual_index:
                g2.write("\n".join([newname2, seq2[read_start:], "+", qual2[read_start:]]) + "\n")
            else:
                g2.write("\n".join([newname2, seq2, "+", qual2]) + "\n")
    return 2 * nseqs


def run_preprocessing(config: PreprocessConfig) -> Tuple[list[str], int]:
    """Start preprocessing."""
    logger.info(f"Start preprocessing of sample {config.sample_name}")

    if config.tmpdir:
        newtmpdir = generate_random_dir(str(config.tmpdir))
    else:
        newtmpdir = generate_random_dir(str(config.output_path))

    # Unzip the fastq.gz files
    if not str(config.read1).endswith("gz"):
        r1file = str(config.read1)
        removerfiles = False
        if config.mode == "paired":
            if not config.read2:
                raise ValueError("Read2 not provided in paired mode")
            r2file = str(config.read2)
    else:
        removerfiles = True
        if config.mode == "paired":
            if not config.read2:
                raise ValueError("Read2 not provided in paired mode")
            r1file = run_unpigz(str(config.read1), newtmpdir, config.num_threads, config.gziptool)
            r2file = run_unpigz(str(config.read2), newtmpdir, config.num_threads, config.gziptool)
        else:
            r1file = run_unpigz(str(config.read1), newtmpdir, config.num_threads, config.gziptool)

    logger.info(f"Writing output files to {config.output_path}")
    # TODO: If fastp was run before, we should skip adapter trimming here.
    if config.adapter_trimming is True:
        output_path = Path(config.output_path)
        r2_input = r2file if config.mode == "paired" else None
        r1file, r2file_out = run_cutadapt(
            r1file, r2_input, output_path, config.sample_name, config.adapter_sequence, config.mode
        )
        if config.mode == "paired":
            if r2file_out is None:
                raise ValueError("Cutadapt returned None for R2 in paired mode")
            r2file = r2file_out

    output_path = Path(config.output_path)
    if config.mode == "single":
        outfilename = str(output_path / f"{config.sample_name}_umis_in_header.fastq")
        nseqs = preprocess_se(r1file, outfilename, config.umi_length, config.spacer_length)
        run_pigz(outfilename, config.num_threads, config.gziptool)
        if removerfiles:
            Path(r1file).unlink()
        Path(newtmpdir).rmdir()
        fastqfiles = [f"{outfilename}.gz"]
    else:
        if config.reverse_index:
            # switch forward and reverse read
            r1filetmp = r1file
            r1file = r2file
            r2file = r1filetmp
            outfile1 = str(output_path / f"{config.sample_name}_R2_umis_in_header.fastq")
            outfile2 = str(output_path / f"{config.sample_name}_R1_umis_in_header.fastq")
        else:
            # r1file=args.read1
            # r2file=args.read2
            outfile1 = str(output_path / f"{config.sample_name}_R1_umis_in_header.fastq")
            outfile2 = str(output_path / f"{config.sample_name}_R2_umis_in_header.fastq")
        nseqs = preprocess_pe(r1file, r2file, outfile1, outfile2, config.umi_length, config.spacer_length, config.dual_index)
        run_pigz(outfile1, config.num_threads, config.gziptool)
        run_pigz(outfile2, config.num_threads, config.gziptool)
        r1path = Path(r1file)
        r2path = Path(r2file)
        if removerfiles is True and r1path.is_file():
            r1path.unlink()
        if removerfiles is True and r2path.is_file():
            r2path.unlink()
        Path(newtmpdir).rmdir()
        fastqfiles = [outfile1 + ".gz", outfile2 + ".gz"]
    logger.info("Finished preprocessing")
    return (fastqfiles, nseqs)
