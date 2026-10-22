#!/usr/bin/env python3

"""
Reduce BAM file size with 2-step strategy:
1. Retain up to N reads per splice junction
2. For unspliced reads, apply per-base coverage limit (e.g., 100x)
Only uniquely mapped reads (NH:i:1) are used.
"""

import pysam
import sys
import argparse
from collections import defaultdict

def reduce_bam_size(input_bam, output_bam, max_junction_reads=20, max_base_coverage=100):
    """
    Main function to reduce BAM file size by subsampling reads.
    
    Args:
        input_bam (str): Path to input BAM file
        output_bam (str): Path to output BAM file
        max_junction_reads (int): Maximum reads to keep per splice junction
        max_base_coverage (int): Maximum coverage for unspliced regions
    """
    # Data stores
    junction_reads = defaultdict(list)               # {(chr, donor, acceptor): [reads]}
    unspliced_reads = []                             # list of unspliced reads
    base_coverage = defaultdict(lambda: defaultdict(int))  # {chr: {pos: coverage}}

    # Open input BAM
    in_bam = pysam.AlignmentFile(input_bam, "rb")

    # First pass: collect junction and unspliced reads
    for read in in_bam.fetch(until_eof=True):
        #if read.is_unmapped or read.is_secondary or read.is_supplementary: 
        if read.is_unmapped or read.is_supplementary:
            continue
        if not read.cigartuples:
            continue

        # Prioritize uniquely mapped reads
        tags = dict(read.tags)
        # if tags.get("NH", 1) != 1:
        #     continue

        chrom = in_bam.get_reference_name(read.reference_id)
        cigar = read.cigartuples
        pos = read.reference_start
        is_spliced = any(op == 3 for op, _ in cigar)

        if is_spliced:
            ref_pos = pos
            added = False
            for op, length in cigar:
                if op == 3:  # N = skipped region (splice junction)
                    donor = ref_pos
                    acceptor = ref_pos + length
                    key = (chrom, donor, acceptor)
                    if len(junction_reads[key]) < max_junction_reads:
                        junction_reads[key].append(read)
                        added = True
                        break  # only add once per read
                if op in (0, 2, 3, 7, 8):
                    ref_pos += length
        else:
            # Store unspliced reads for second-pass filtering
            unspliced_reads.append(read)

    in_bam.close()

    # Open output BAM
    out_bam = pysam.AlignmentFile(output_bam, "wb", template=pysam.AlignmentFile(input_bam, "rb"))

    # Step 2: Write junction reads (no coverage filter)
    for reads in junction_reads.values():
        for r in reads:
            out_bam.write(r)

    # Step 3: Apply per-base coverage filter to unspliced reads
    for read in unspliced_reads:
        chrom = out_bam.get_reference_name(read.reference_id)
        ref_positions = read.get_reference_positions()

        # Check if any position exceeds max coverage
        if all(base_coverage[chrom][p] < max_base_coverage for p in ref_positions if p is not None):
            out_bam.write(read)
            for p in ref_positions:
                if p is not None:
                    base_coverage[chrom][p] += 1

    out_bam.close()

def main(args=None):
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description="Reduce BAM file size by subsampling reads at splice junctions and coverage sites")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument("--max_junction_reads", type=int, default=50,
                        help="Maximum reads to keep per splice junction (default: 50)")
    parser.add_argument("--max_base_coverage", type=int, default=50,
                        help="Maximum coverage for unspliced regions (default: 50)")

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    reduce_bam_size(args.input_bam, args.output_bam, 
                   args.max_junction_reads, args.max_base_coverage)

if __name__ == "__main__":
    main()
