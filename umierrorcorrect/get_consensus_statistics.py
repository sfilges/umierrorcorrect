#!/usr/bin/env python3
import logging
import random
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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


def plot_downsampling(results_tot, fsize, plot_filename):
    x = []
    y = []
    for r in results_tot[0]:
        h = results_tot[0][r]
        x.append(h.total_reads[int(fsize)])
        y.append(h.umis[int(fsize)])
    plt.plot(x, y, "o-")
    plt.xlabel("Depth")
    plt.ylabel("Number of UMI families")
    plt.title("Downsampling plot")
    plt.box(False)
    plt.xlim(0, max(x) + 40000)
    plt.ylim(0, max(y) + 1000)
    plt.savefig(plot_filename)


def save_downsampled_table(all_results, tot_results, out_filename):
    with Path(out_filename).open("w") as g:
        for r in tot_results[0]:
            text = tot_results[0][r].write_stats()
            lines = text.split("\n")
            for line in lines:
                g.write("downsampled" + str(r) + "\t" + line + "\n")
        for region in all_results:
            for r in region:
                text = region[r].write_stats()
                lines = text.split("\n")
                for line in lines:
                    g.write("downsampled" + str(r) + "\t" + line + "\n")


def downsample_reads_per_region(hist, _fraction, fsizes, onlyNamed=True):
    all_results = []
    downsample_rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for h in hist:
        run_analysis = True
        if onlyNamed and h.name == "":
            run_analysis = False
        if run_analysis:
            num_families = len(h.hist)
            tmpnames = np.array(range(0, num_families))
            singnames = list(range(num_families, num_families + h.singletons))
            times = np.array(h.hist)
            reads = np.repeat(tmpnames, times, axis=0)  # expand to one entry per read
            reads = list(reads) + singnames
            results = {}
            for r in downsample_rates:
                ds_reads = random.sample(list(reads), round(r * len(reads)))  # noqa: S311 - downsample
                new_hist = Counter(ds_reads).values()  # collapse to one entry per UMI family
                new_hist = sorted(new_hist, reverse=True)  # sort
                new_singletons = list(new_hist).count(1)  # count singletons in new
                new_stat = region_cons_stat(h.regionid, h.pos, h.name, new_singletons, h.fsizes)
                new_stat.add_histogram(new_hist, fsizes)
                results[r] = new_stat
            all_results.append(results)
    return all_results


# def write_report():
#    from markdown2 import Markdown
#    md=Markdown()
#    print(md.convert("#Report for sample {}".format("samplename")))
#
#
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
