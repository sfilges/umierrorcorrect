#!/usr/bin/env python3
"""Consensus statistics functionality for UMI error correction.

This module provides functions for calculating and reporting statistics
from consensus BAM files and histogram data.
"""
import logging
from collections import Counter
from pathlib import Path

import pysam


class region_cons_stat:
    def __init__(self, regionid, pos, name, singletons, fsizes):
        self.regionid = regionid
        self.pos = pos
        self.name = name
        self.singletons = singletons
        self.hist = []
        self.total_reads = {}
        self.umis = {}
        for fsize in fsizes:
            self.total_reads[fsize] = 0
            self.umis[fsize] = 0
            self.total_reads[0] = self.singletons
            self.umis[0] = self.singletons
        self.total_reads[1] = self.singletons
        self.umis[1] = self.singletons
        self.fsizes = fsizes

    def add_histogram(self, hist, fsizes):
        self.total_reads[0] += sum(hist)
        self.umis[0] += sum(hist)
        for fsize in fsizes:
            tmp = [x for x in hist if x >= fsize]
            if fsize == 1:
                self.total_reads[fsize] += sum(tmp)
                self.umis[fsize] += len(tmp)
            else:
                self.total_reads[fsize] += sum(tmp)
                self.umis[fsize] += len(tmp)
        self.hist.extend(hist)

    def write_stats(self):
        lines = []
        r0 = self.total_reads[0]
        u0 = self.umis[0]
        line = "\t".join([str(self.regionid), self.pos, self.name, "0", "1.0", str(r0), str(u0)])
        lines.append(line)
        for fsize in self.fsizes:
            fraction = 0 if r0 == 0 else self.total_reads[fsize] / r0
            line = "\t".join(
                [
                    str(self.regionid),
                    self.pos,
                    self.name,
                    str(fsize),
                    str(1.0 * fraction),
                    str(self.total_reads[fsize]),
                    str(self.umis[fsize]),
                ]
            )
            lines.append(line)
        return "\n".join(lines)


def get_stat(consensus_filename, stat_filename):
    with Path(stat_filename).open() as f:
        regions = []
        for line in f:
            line = line.rstrip()
            regionid, pos, name, cons, singles, *rest = line.split("\t")
            singles = int(singles.split(": ")[-1])
            regionid = str(regionid)
            regions.append((regionid, pos, singles, name))

    # print(regions)
    hist = {}
    with pysam.AlignmentFile(consensus_filename, "rb") as f:
        reads = f.fetch()
        for read in reads:
            idx = read.qname
            if idx.startswith("Consensus_read"):
                parts = idx.split("_")
                regionid = str(parts[2])
                # regionid='_'.join(parts[2:-2])
                if parts[-1].startswith("Count") or parts[-1] == "a":
                    count = int(idx.split("=")[-1])
                    if regionid not in hist:
                        hist[regionid] = []
                    hist[regionid].append(count)
    # print(hist)
    fsizes = [1, 2, 3, 4, 5, 7, 10, 20, 30]
    regionstats = []
    for regionid, pos, singletons, name in regions:
        if "-" in regionid:
            a, b, *rest = regionid.split("-")
            from_tag = False
            try:
                int(b)
            except ValueError:  # not an int
                if "_" in b:
                    name = a
                    a = 0
                    b = int(b.split("_")[-1])
                    from_tag = True

            stat = region_cons_stat(regionid, pos, name, singletons, fsizes)
            for i in range(int(a), int(b) + 1):
                if not from_tag:
                    if str(i) in hist:
                        stat.add_histogram(hist[str(i)], fsizes)
                else:
                    if name + "_" + str(i) in hist:
                        stat.add_histogram(hist[name + "_" + str(i)], fsizes)
            # print(stat.write_stats())
            regionstats.append(stat)
        else:
            stat = region_cons_stat(regionid, pos, name, singletons, fsizes)
            if regionid in hist:
                stat.add_histogram(hist[regionid], fsizes)
            regionstats.append(stat)
    return regionstats


def calculate_target_coverage(stats, fsizes):
    reads_all = {}
    reads_target = {}
    fsizes.insert(0, 0)
    for fsize in fsizes:
        reads_all[fsize] = 0
        reads_target[fsize] = 0
    for region in stats:
        for fsize in fsizes:
            reads_all[fsize] += region.umis[fsize]
            if region.name not in "":
                reads_target[fsize] += region.umis[fsize]
    lines = []
    for fsize in fsizes:
        if reads_all[fsize] > 0:
            lines.append(
                f"{fsize}\t{reads_target[fsize]}\t{reads_all[fsize]}\t{1.0 * (reads_target[fsize] / reads_all[fsize])}"
            )
        else:
            lines.append(f"{fsize}\t{reads_target[fsize]}\t{reads_all[fsize]}\t{0}")

    return "\n".join(lines)


def get_overall_statistics(hist, fsizes):
    histall = region_cons_stat("All", "all_regions", "", 0, fsizes)
    fsizesnew = fsizes.copy()
    histall.fsizes = fsizes
    fsizesnew.insert(0, 0)
    # print(fsizesnew)
    for fsize in fsizesnew:
        histall.total_reads[fsize] = 0
        histall.umis[fsize] = 0

    for region in hist:
        for fsize in fsizesnew:
            histall.total_reads[fsize] += region.total_reads[fsize]
            histall.umis[fsize] += region.umis[fsize]
    # print(histall.write_stats())
    return histall


def get_percent_mapped_reads(num_fastq_reads, bamfile):
    """Get the number of mapped reads from the .bai file"""
    with pysam.AlignmentFile(bamfile, "rb") as f:
        stats = f.get_index_statistics()
        num_mapped = 0
        for s in stats:
            num_mapped += s.mapped
    ratio = (num_mapped / num_fastq_reads) * 1.0
    return (num_mapped, ratio)


def run_get_consensus_statistics(output_path, consensus_filename, stat_filename, output_raw, samplename):
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
    histall = get_overall_statistics(hist, fsizes)
    if not consensus_filename:
        consensus_filename = str(list(out_path.glob("*_consensus_reads.bam"))[0])
        # print(consensus_filename)
    if not samplename:
        samplename = Path(stat_filename).stem
    outfilename = out_path / f"{samplename}_summary_statistics.txt"
    logging.info(f"Writing consensus statistics to {outfilename}")
    with outfilename.open("w") as g:
        g.write(histall.write_stats() + "\n")
        for stat in hist:
            g.write(stat.write_stats() + "\n")
    outfilename = out_path / f"{samplename}_target_coverage.txt"
    with outfilename.open("w") as g:
        g.write(calculate_target_coverage(hist, fsizes))
    if output_raw:
        largehist = []
        outfilename = out_path / f"{samplename}_consensus_group_counts.txt"
        for h in hist:
            largehist = largehist + h.hist
            largehist = largehist + [1] * h.singletons
        hist_counts = Counter(largehist)
        # print(hist_counts)
        with outfilename.open("w") as g:
            for size in sorted(hist_counts):
                g.write(str(size) + "\t" + str(hist_counts[size]) + "\n")

    # print(hist)
    # plot_histogram(hist,output_path+'/histogram.png')
    logging.info("Finished consensus statistics")
    # write_report()


def main(output_path, consensus_filename, stat_filename, output_raw, samplename):
    run_get_consensus_statistics(output_path, consensus_filename, stat_filename, output_raw, samplename)
