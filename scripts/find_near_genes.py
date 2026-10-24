#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import sys
import os

def clean_chr(chr_val):
    """ chr01 -> chr1, 12 -> chr12"""
    c = str(chr_val).strip()
    if c.startswith("chr"):
        c = c[3:]
    try:
        return f"chr{int(c)}"
    except ValueError:
        return f"chr{c}"

def parse_args():
    parser = argparse.ArgumentParser(description="Find genes within a specified window near hotspot SNPs")
    parser.add_argument("--hotspot_file", type=str, required=True,
                        help="Input hotspot SNP file (e.g., meta_00400_points_hotspot.txt)")
    parser.add_argument("--gene_loc", type=str, required=True,
                        help="Gene location file (e.g., DM8.1_gene.location.txt)")
    parser.add_argument("--dist", type=float, default=2.84,
                        help="Upstream and downstream extension distance from SNP, in kb (default: 2.84, i.e., ±2840 bp)")
    parser.add_argument("--output", type=str, default="snps_with_nearby_genes.txt",
                        help="Output result filename (default: snps_with_nearby_genes.txt)")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Convert distance to bp (e.g., 2.84 kb -> 2840 bp)
    dist_bp = int(args.dist * 1000)
    
    print(f"==================================================")
    print(f" Hotspot SNP file: {args.hotspot_file}")
    print(f" Gene location file: {args.gene_loc}")
    print(f" Search window: ±{args.dist} kb (±{dist_bp} bp)")
    print(f"==================================================\n")

    # 1. Read gene location file (no header)
    if not os.path.exists(args.gene_loc):
        print(f"[Error] Gene location file does not exist: {args.gene_loc}", file=sys.stderr)
        sys.exit(1)

    gene_df = pd.read_csv(args.gene_loc, sep=r'\s+', header=None,
                          names=['Gene_ID', 'Chr', 'Gene_Start', 'Gene_End', 'Strand'])
    
    # Clean and convert gene coordinate data types
    gene_df['Gene_Start'] = pd.to_numeric(gene_df['Gene_Start'], errors='coerce')
    gene_df['Gene_End'] = pd.to_numeric(gene_df['Gene_End'], errors='coerce')
    gene_df['clean_chr'] = gene_df['Chr'].apply(clean_chr)
    gene_df = gene_df.dropna(subset=['Gene_Start', 'Gene_End']).copy()
    gene_df['Gene_Start'] = gene_df['Gene_Start'].astype(int)
    gene_df['Gene_End'] = gene_df['Gene_End'].astype(int)

    # 2. Read SNP file
    if not os.path.exists(args.hotspot_file):
        print(f"[Error] SNP file does not exist: {args.hotspot_file}", file=sys.stderr)
        sys.exit(1)

    snp_df = pd.read_csv(args.hotspot_file, sep=r'\s+')
    
    # Match chromosome and BP columns
    col_map = {col.upper(): col for col in snp_df.columns}
    chr_col = col_map.get('CHR')
    bp_col = col_map.get('SNP_BP') or col_map.get('BP') or col_map.get('POS')

    if not (chr_col and bp_col):
        print(f"[Error] Cannot identify Chr or SNP_BP column in the SNP file!", file=sys.stderr)
        sys.exit(1)

    snp_df['clean_chr'] = snp_df[chr_col].apply(clean_chr)
    snp_df[bp_col] = pd.to_numeric(snp_df[bp_col], errors='coerce')
    snp_df = snp_df.dropna(subset=[bp_col]).copy()
    snp_df[bp_col] = snp_df[bp_col].astype(int)

    # 3. Match genes within 2.84kb for each SNP
    results = []

    for _, snp_row in snp_df.iterrows():
        s_chr = snp_row['clean_chr']
        s_bp = snp_row[bp_col]
        
        # Define search window [s_bp - dist_bp, s_bp + dist_bp]
        win_start = max(0, s_bp - dist_bp)
        win_end = s_bp + dist_bp

        # Filter genes on the same chromosome that overlap with the window
        # Overlap condition: Gene_Start <= win_end AND Gene_End >= win_start
        matched_genes = gene_df[
            (gene_df['clean_chr'] == s_chr) & 
            (gene_df['Gene_Start'] <= win_end) & 
            (gene_df['Gene_End'] >= win_start)
        ].copy()

        if not matched_genes.empty:
            for _, g_row in matched_genes.iterrows():
                # Calculate the nearest physical distance from SNP to gene boundaries (0 if SNP is inside the gene)
                if s_bp < g_row['Gene_Start']:
                    dist_to_gene = g_row['Gene_Start'] - s_bp
                elif s_bp > g_row['Gene_End']:
                    dist_to_gene = s_bp - g_row['Gene_End']
                else:
                    dist_to_gene = 0

                res_dict = snp_row.to_dict()
                # Remove helper column
                res_dict.pop('clean_chr', None)
                
                # Append gene information
                res_dict['Nearby_Gene_ID'] = g_row['Gene_ID']
                res_dict['Gene_Chr'] = g_row['Chr']
                res_dict['Gene_Start'] = g_row['Gene_Start']
                res_dict['Gene_End'] = g_row['Gene_End']
                res_dict['Gene_Strand'] = g_row['Strand']
                res_dict['Dist_to_Gene_bp'] = dist_to_gene
                
                results.append(res_dict)
        else:
            # If no matching gene is found nearby, keep the SNP record and mark as empty
            res_dict = snp_row.to_dict()
            res_dict.pop('clean_chr', None)
            res_dict['Nearby_Gene_ID'] = "None"
            res_dict['Gene_Chr'] = "-"
            res_dict['Gene_Start'] = "-"
            res_dict['Gene_End'] = "-"
            res_dict['Gene_Strand'] = "-"
            res_dict['Dist_to_Gene_bp'] = "-"
            results.append(res_dict)

    # 4. Export results
    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"[Done] Analysis finished. Found {len(out_df[out_df['Nearby_Gene_ID'] != 'None'])} matched associations!")
    print(f"Results exported to: {args.output}\n")
    
    # Preview the first few rows in the terminal
    print("=== Result Preview (First 5 rows) ===")
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    print(out_df.head().to_string(index=False))

if __name__ == "__main__":
    main()