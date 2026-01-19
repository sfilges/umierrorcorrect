#!/usr/bin/env python3
import logging
from pathlib import Path

from umierrorcorrect.get_consensus_statistics import (
    downsample_reads_per_region,
    get_stat,
    plot_downsampling,
    region_cons_stat,
    save_downsampled_table,
)


def run_downsampling(output_path, consensus_filename, stat_filename, fsize, samplename=None):
    logging.info("Getting consensus statistics")
    out_path = Path(output_path)
    if not consensus_filename:
        consensus_filename = str(list(out_path.glob("*_consensus_reads.bam"))[0])
    if not samplename:
        samplename = Path(consensus_filename).name.replace("_consensus_reads.bam", "")
    if not stat_filename:
        stat_filename = str(out_path / f"{samplename}.hist")
    hist = get_stat(consensus_filename, stat_filename)
    fsizes = [1, 2, 3, 4, 5, 7, 10, 20, 30]
    tot_results = region_cons_stat("All", "all_regions", "", 0, fsizes)
    for h in hist:
        tot_results.hist = tot_results.hist + h.hist
        tot_results.singletons += h.singletons

    downsample_rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    tot = downsample_reads_per_region([tot_results], downsample_rates, fsizes, False)
    all_results = downsample_reads_per_region(hist, downsample_rates, fsizes, True)
    filename = str(out_path / f"{samplename}_downsampled_coverage.txt")
    save_downsampled_table(all_results, tot, filename)
    filename = str(out_path / f"{samplename}_downsampled_plot.png")
    plot_downsampling(tot, fsize, filename)
