#!/usr/bin/env python3
import pysam


def filter_bam(infilename, outfilename, consensus_cutoff):
    consensus_cutoff = int(consensus_cutoff)
    with pysam.AlignmentFile(infilename, "rb") as f, pysam.AlignmentFile(outfilename, "wb", template=f) as g:
        reads = f.fetch()
        for read in reads:
            size = int(read.qname.rsplit("=", 1)[-1])
            if size >= consensus_cutoff:
                g.write(read)
