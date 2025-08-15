#!/usr/bin/env python3
"""
Preprocess lncRNA data to extract proper sample names and clean the format.
"""

import pandas as pd
import re

def preprocess_lncrna_data():
    """Preprocess the lncRNA data file."""
    print("Preprocessing lncRNA data...")
    
    # Read the file line by line to find the header
    with open('../data/Apul_lncRNA_counts_filtered.txt', 'r') as f:
        lines = f.readlines()
    
    # Find the line with column headers (contains Geneid, Chr, Start, etc.)
    header_line_idx = None
    for i, line in enumerate(lines):
        if line.startswith('Geneid\tChr\tStart\tEnd\tStrand\tLength'):
            header_line_idx = i
            break
    
    if header_line_idx is None:
        print("Error: Could not find header line")
        return None
    
    print(f"Found header at line {header_line_idx + 1}")
    
    # Read the data starting from the header line
    lncrna_data = pd.read_csv('../data/Apul_lncRNA_counts_filtered.txt', 
                             sep='\t', 
                             skiprows=header_line_idx)
    
    print(f"Loaded {lncrna_data.shape[0]} lncRNAs with {lncrna_data.shape[1]} columns")
    
    # Find count columns (those with ACR- pattern)
    count_cols = [col for col in lncrna_data.columns if 'ACR-' in col]
    print(f"Found {len(count_cols)} count columns")
    
    # Clean column names - extract just the ACR-XXX-TPX part
    clean_cols = {}
    for col in count_cols:
        if 'ACR-' in col:
            # Extract ACR-XXX-TPX pattern
            match = re.search(r'ACR-\d+-TP\d+', col)
            if match:
                clean_name = match.group()
                clean_cols[col] = clean_name
                print(f"  {col} -> {clean_name}")
    
    # Rename columns
    lncrna_data = lncrna_data.rename(columns=clean_cols)
    
    # Keep only essential columns
    essential_cols = ['Geneid'] + list(clean_cols.values())
    lncrna_data_clean = lncrna_data[essential_cols].copy()
    
    # Set Geneid as index
    lncrna_data_clean.set_index('Geneid', inplace=True)
    
    # Remove rows with all zero counts
    lncrna_data_clean = lncrna_data_clean.loc[(lncrna_data_clean != 0).any(axis=1)]
    
    print(f"Final dataset: {lncrna_data_clean.shape[0]} lncRNAs, {lncrna_data_clean.shape[1]} samples")
    
    # Save cleaned data
    output_file = '../data/lncrna_counts_cleaned.csv'
    lncrna_data_clean.to_csv(output_file)
    print(f"Cleaned data saved to {output_file}")
    
    return lncrna_data_clean

if __name__ == "__main__":
    preprocess_lncrna_data()
