#!/usr/bin/env python3
"""
filter_diamond.py

Filters a diamond TSV report for viral contigs and applies a length threshold
based on the 95th percentile (with a minimum of 300 bp and a maximum cap of 900 bp).
The output CSV filename is set to 'diamond_parsing.csv'.
Additionally, a text file 'length_threshold.txt' records the threshold value used.

Usage:
    python filter_diamond.py --input diamond.tsv
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description='Filter diamond report for Viruses and apply length threshold.'
    )
    parser.add_argument(
        '--input', '-i', required=True,
        help='Path to the diamond TSV file'
    )
    parser.add_argument(
        '--percentile', '-p', type=float, default=95.0,
        help='Percentile for length threshold (default: 95)'
    )
    parser.add_argument(
        '--min_length', '-m', type=int, default=300,
        help='Minimum enforced threshold if percentile is lower (default: 300)'
    )
    parser.add_argument(
        '--max_length', '-x', type=int, default=900,
        help='Maximum enforced threshold if percentile is higher (default: 900)'
    )
    args = parser.parse_args()

    diamond_file = args.input
    if not os.path.isfile(diamond_file):
        print(f"Error: file '{diamond_file}' not found.")
        sys.exit(1)

    # 1. Read diamond TSV file into DataFrame
    df = pd.read_csv(diamond_file, sep='\t', low_memory=False)

    # 2. Ensure required columns exist and enforce correct data types
    required_cols = ['bitscore', 'qseqid']
    for col in required_cols:
        if col not in df.columns:
            print(f"Error: '{col}' column not found in input.")
            sys.exit(1)

    # 2a. Convert 'evalue' to numeric
    # try:
    #     df['evalue'] = pd.to_numeric(df['evalue'], errors='raise')
    # except Exception:
    #     print("Error: Failed to convert 'evalue' to numeric. Please check input format.")
    #     sys.exit(1)

    # 2b. Convert 'bitscore' to numeric
    try:
        df['bitscore'] = pd.to_numeric(df['bitscore'], errors='raise')
    except Exception:
        print("Error: Failed to convert 'bitscore' to numeric. Please check input format.")
        sys.exit(1)

    # 2c. Convert 'qseqid' to string
    try:
        df['qseqid'] = df['qseqid'].astype(str)
    except Exception:
        print("Error: Failed to convert 'qseqid' to string. Please check input format.")
        sys.exit(1)

    # 3. Filter alignments by e-value threshold
    # df = df[df['evalue'] < 1e-5]

    # 4. For each qseqid, keep only the row with the highest bitscore
    df = (
        df
        .sort_values('bitscore', ascending=False)
        .drop_duplicates(subset=['qseqid'], keep='first')
    )

    # 5. Filter for kingdom 'Viruses'
    if 'sskingdoms' not in df.columns:
        print("Error: 'sskingdoms' column not found.")
        sys.exit(1)
    virus_df = df[df['sskingdoms'] == 'Viruses']
    virus_df.to_csv('diamond_virus.csv', index=False)

    # 6. Check for 'length' column and calculate percentile threshold
    if 'length' not in virus_df.columns:
        print("Error: 'length' column not found in diamond_virus.csv.")
        sys.exit(1)
    lengths = virus_df['length'].dropna()
    if lengths.empty:
        print("Error: No length values available for percentile calculation.")
        sys.exit(1)

    threshold = np.percentile(lengths, args.percentile)

    # 7. Enforce minimum and maximum threshold boundaries
    if threshold < args.min_length:
        threshold = args.min_length
        if not (lengths > threshold).any():
            print(f"No sequences longer than {args.min_length}")
            sys.exit(0)
    elif threshold > args.max_length:
        threshold = args.max_length
        if not (lengths > threshold).any():
            print(f"No sequences longer than {args.max_length}")
            sys.exit(0)

    # 8. Filter sequences that exceed the final threshold
    filtered_df = virus_df[virus_df['length'] > threshold]

    # 9. Save filtered results to CSV with fixed filename
    output_file = 'diamond_parsing.csv'
    filtered_df.to_csv(output_file, index=False)
    print(f"Filtered sequences saved to '{output_file}' with threshold = {threshold}")

    # 10. Record the threshold value to a text file
    try:
        with open('length_threshold.txt', 'w') as tf:
            tf.write(str(threshold))
    except Exception:
        print("Warning: Could not write threshold to 'length_threshold.txt'.")

if __name__ == '__main__':
    main()
