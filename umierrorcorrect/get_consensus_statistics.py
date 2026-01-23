#!/usr/bin/env python3
"""Consensus statistics functionality for UMI error correction.

This module provides functions for calculating and reporting statistics
from consensus BAM files and histogram data.
"""

from collections import Counter
from pathlib import Path

import pysam

from umierrorcorrect.core.constants import DEFAULT_FAMILY_SIZES, HISTOGRAM_SUFFIX
from umierrorcorrect.core.logging_config import get_logger

logger = get_logger(__name__)


class RegionConsensusStats:
    def __init__(self, regionid, pos, name, singletons, fsizes):
        self.regionid = regionid
        self.pos = pos
        self.name = name
        self.singletons = singletons
        self.family_sizes = []
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

    def add_family_sizes(self, family_sizes, fsizes):
        # Threshold 0 represents raw read statistics (no UMI deduplication)
        # Total reads and UMIs are both equal to the sum of all reads in the families
        self.total_reads[0] += sum(family_sizes)
        self.umis[0] += sum(family_sizes)

        # Thresholds >= 1 represent unique molecular statistics
        for fsize in fsizes:
            tmp = [x for x in family_sizes if x >= fsize]
            self.total_reads[fsize] += sum(tmp)
            self.umis[fsize] += len(tmp)
        self.family_sizes.extend(family_sizes)

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


def get_stat(consensus_filename: Path, stat_filename: Path) -> list:
    """
    Get consensus statistics per region from consensus BAM file and stat file.

    Parameters
    ----------
    consensus_filename : Path
        Path to consensus BAM file.
    stat_filename : Path
        Path to stat file.

    Returns
    -------
    list
        List of region statistics.
    """

    # ---------------------------------------------------------------
    # Get regions from stat file
    # ---------------------------------------------------------------
    # TODO: Should this include only regions with names if any are present, i.e. if a bed file was used? On the other hand
    # non-annotated regions contain off-target reads which may be interesting to retain for downstream analysis.
    with Path(stat_filename).open() as f:
        regions = []
        for line in f:
            line = line.rstrip()
            regionid, pos, name, cons, singles, *rest = line.split("\t")
            singles = int(singles.split(": ")[-1])
            regionid = str(regionid)
            regions.append((regionid, pos, singles, name))

    # ---------------------------------------------------------------
    # Read consensus BAM file and get consensus read histograms
    # ---------------------------------------------------------------

    # Read names have the structure: type_read_regionid_umi_Count=number
    # singletons: Singleton_read_7_GTCGAAACTAGA_Count=1
    # consensus: Consensus_read_7_TGATAAAATAAG_a_Count=8
    # What do the a, b, c, etc. tags before counts mean? Likely subclusters from UMI clustering.
    family_sizes_by_region: dict[str, list[int]] = {}
    with pysam.AlignmentFile(consensus_filename, "rb") as f:
        reads = f.fetch()
        for read in reads:
            idx = read.qname  # type: ignore
            if idx.startswith("Consensus_read"):
                # ['Consensus', 'read', '7', 'TGATAAAATAAG', 'a', 'Count=8']
                parts = idx.split("_")
                regionid = str(parts[2])
                if parts[-1].startswith("Count") or parts[-1] == "a":
                    count = int(idx.split("=")[-1])
                    if regionid not in family_sizes_by_region:
                        family_sizes_by_region[regionid] = []
                    # Dictionary with regionid as key and list of counts as value
                    # family_sizes_by_region = {regionid: [count, count, count, ...]}
                    # len(family_sizes_by_region[regionid]) = number of consensus reads per region
                    family_sizes_by_region[regionid].append(count)

    fsizes = list(DEFAULT_FAMILY_SIZES)[1:]  # Exclude 0, which is handled separately
    regionstats = []
    for regionid, pos, singletons, name in regions:
        if "-" in regionid:
            a_str, b_str, *rest = regionid.split("-")
            from_tag = False
            start_idx = 0
            end_idx = 0
            try:
                start_idx = int(a_str)
                end_idx = int(b_str)
            except ValueError:  # not an int
                if "_" in b_str:
                    name = a_str
                    start_idx = 0
                    end_idx = int(b_str.split("_")[-1])
                    from_tag = True

            stat = RegionConsensusStats(regionid, pos, name, singletons, fsizes)
            for i in range(start_idx, end_idx + 1):
                if not from_tag:
                    if str(i) in family_sizes_by_region:
                        stat.add_family_sizes(family_sizes_by_region[str(i)], fsizes)
                else:
                    if name + "_" + str(i) in family_sizes_by_region:
                        stat.add_family_sizes(family_sizes_by_region[name + "_" + str(i)], fsizes)
            regionstats.append(stat)
        else:
            stat = RegionConsensusStats(regionid, pos, name, singletons, fsizes)
            if regionid in family_sizes_by_region:
                stat.add_family_sizes(family_sizes_by_region[regionid], fsizes)
            regionstats.append(stat)
    return regionstats


def calculate_target_coverage(stats, fsizes=None):
    if fsizes is None:
        fsizes = stats[0].fsizes if stats else list(DEFAULT_FAMILY_SIZES)[1:]

    reads_all = {}
    reads_target = {}
    fsizes_calc = fsizes.copy()
    fsizes_calc.insert(0, 0)
    for fsize in fsizes_calc:
        reads_all[fsize] = 0
        reads_target[fsize] = 0
    for region in stats:
        for fsize in fsizes_calc:
            reads_all[fsize] += region.umis[fsize]
            if region.name not in "":
                reads_target[fsize] += region.umis[fsize]
    lines = []
    for fsize in fsizes_calc:
        if reads_all[fsize] > 0:
            lines.append(
                f"{fsize}\t{reads_target[fsize]}\t{reads_all[fsize]}\t{1.0 * (reads_target[fsize] / reads_all[fsize])}"
            )
        else:
            lines.append(f"{fsize}\t{reads_target[fsize]}\t{reads_all[fsize]}\t{0}")

    return "\n".join(lines)


def get_overall_statistics(hist, fsizes=None):
    if fsizes is None:
        fsizes = hist[0].fsizes if hist else list(DEFAULT_FAMILY_SIZES)[1:]

    histall = RegionConsensusStats("All", "all_regions", "", 0, fsizes)
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
    logger.info("Getting consensus statistics")
    out_path = Path(output_path)
    if not consensus_filename:
        consensus_filename = str(list(out_path.glob("*_consensus_reads.bam"))[0])
    if not samplename:
        samplename = Path(consensus_filename).name.replace("_consensus_reads.bam", "")
    if not stat_filename:
        stat_filename = str(out_path / f"{samplename}{HISTOGRAM_SUFFIX}")
    hist = get_stat(consensus_filename, stat_filename)
    fsizes = list(DEFAULT_FAMILY_SIZES)[1:]  # Exclude 0, which is handled separately
    histall = get_overall_statistics(hist, fsizes)
    if not consensus_filename:
        consensus_filename = str(list(out_path.glob("*_consensus_reads.bam"))[0])
        # print(consensus_filename)
    if not samplename:
        samplename = Path(stat_filename).stem
    outfilename = out_path / f"{samplename}_summary_statistics.txt"
    logger.info(f"Writing consensus statistics to {outfilename}")
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

    logger.info("Finished consensus statistics")


def main(output_path, consensus_filename, stat_filename, output_raw, samplename):
    run_get_consensus_statistics(output_path, consensus_filename, stat_filename, output_raw, samplename)
