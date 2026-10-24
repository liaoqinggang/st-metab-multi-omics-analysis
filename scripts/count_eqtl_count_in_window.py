#!/usr/bin/env python3
"""
Script Purpose: Based on genome chromosome lengths, divide each chromosome into fixed-size windows 
(e.g., 1000kb) and count the number of SNP loci in each window from cis and trans files whose 
absolute correlation value exceeds a given threshold.

Usage:
python3 count_snp_in_windows.py <genome_len_file> <cis_file> <trans_file> [threshold] [window_size]

Parameters:
    genome_len_file : Genome length file (e.g., DM8.1_genome.len.txt)
    cis_file        : Cis SNP file (e.g., cis_snp_r_00400_pathway.txt)
    trans_file      : Trans SNP file (e.g., trans_snp_r_00400.txt)
    threshold       : Absolute correlation threshold, default 0.5
    window_size     : Window size in bp, default 1000000 (1000kb)

Output: Tab-separated 5-column result
    Chromosome ID   Window start position   Window end position   Cis count   Trans count
"""

import sys
import os
from collections import defaultdict

def load_genome_lengths(genome_file):
    """Load genome length file, keeping only numbered chromosomes (chr01~chr12, etc.), filtering out contigs/scaffolds."""
    chrom_lengths = {}
    with open(genome_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                length = int(parts[1])
            except ValueError:
                continue
            # Keep only chromosomes starting with 'chr' followed by digits, filtering contigs/scaffolds
            if chrom.startswith('chr') and chrom[3:].isdigit():
                chrom_lengths[chrom] = length
    return chrom_lengths

def load_snp_positions(snp_file, threshold):
    """Read SNP file, parse chromosome and position, and filter loci with absolute correlation greater than threshold."""
    snp_positions = defaultdict(list)  # Use defaultdict to auto-initialize empty lists and avoid KeyError
    with open(snp_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            snp_id = parts[0]
            try:
                corr_value = float(parts[5])
            except ValueError:
                continue
            # Filter loci where absolute correlation value is greater than threshold
            if abs(corr_value) <= threshold:
                continue
            # Parse chromosome and position from SNP ID, e.g., chr01_71886511
            if '_' not in snp_id:
                continue
            chrom_pos = snp_id.rsplit('_', 1)
            if len(chrom_pos) != 2:
                continue
            chrom, pos_str = chrom_pos
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            snp_positions[chrom].append(pos)  # defaultdict automatically handles the first appearance
    return dict(snp_positions)  # Convert back to standard dict for subsequent compatibility

def count_snps_in_windows(chrom_lengths, cis_positions, trans_positions, window_size):
    """Count the number of cis and trans loci within each window across chromosomes."""
    results = []
    for chrom in sorted(chrom_lengths.keys()):
        length = chrom_lengths[chrom]
        cis_list = sorted(cis_positions.get(chrom, []))
        trans_list = sorted(trans_positions.get(chrom, []))
        
        # Sliding window
        window_start = 0
        while window_start < length:
            window_end = min(window_start + window_size, length)
            
            # Count cis loci
            cis_count = sum(1 for pos in cis_list if window_start <= pos < window_end)
            # Count trans loci
            trans_count = sum(1 for pos in trans_list if window_start <= pos < window_end)
            
            results.append((chrom, window_start, window_end, cis_count, trans_count))
            window_start += window_size
            
    return results

def main():
    if len(sys.argv) < 4:
        print("Usage: python3 count_snp_in_windows.py <genome_len_file> <cis_file> <trans_file> [threshold] [window_size]")
        print("Example: python3 count_snp_in_windows.py DM8.1_genome.len.txt cis_snp_r_00400_pathway.txt trans_snp_r_00400.txt 0.5 20000")
        sys.exit(1)
        
    genome_file = sys.argv[1]
    cis_file = sys.argv[2]
    trans_file = sys.argv[3]
    threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    window_size = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
        
    # Check if files exist
    for f in [genome_file, cis_file, trans_file]:
        if not os.path.exists(f):
            print(f"Error: File does not exist - {f}")
            sys.exit(1)
        
    # Load data
    print("Loading genome length file...", file=sys.stderr)
    chrom_lengths = load_genome_lengths(genome_file)
    print(f"  Loaded {len(chrom_lengths)} chromosomes", file=sys.stderr)
        
    print("Loading cis SNP file...", file=sys.stderr)
    cis_positions = load_snp_positions(cis_file, threshold)
    cis_total = sum(len(v) for v in cis_positions.