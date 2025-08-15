#!/usr/bin/env python3
"""
Preprocess all data files to ensure consistent sample naming and clean format.
"""

import pandas as pd
import re

def preprocess_mirna_data():
    """Preprocess miRNA data."""
    print("Preprocessing miRNA data...")
    
    mirna_data = pd.read_csv('../data/Apul_miRNA_counts_formatted.txt', sep='\t')
    
    # Clean column names - extract just the ACR-XXX-TPX part
    clean_cols = {}
    for col in mirna_data.columns:
        if 'ACR-' in col and 'TP' in col:
            # Extract ACR-XXX-TPX pattern from various formats
            match = re.search(r'ACR-\d+_TP\d+', col)
            if match:
                clean_name = match.group().replace('_', '-')
                clean_cols[col] = clean_name
                print(f"  {col} -> {clean_name}")
    
    # Rename columns
    mirna_data = mirna_data.rename(columns=clean_cols)
    
    # Keep only essential columns
    essential_cols = ['Name'] + list(clean_cols.values())
    mirna_data_clean = mirna_data[essential_cols].copy()
    
    # Set Name as index
    mirna_data_clean.set_index('Name', inplace=True)
    
    # Remove rows with all zero counts
    mirna_data_clean = mirna_data_clean.loc[(mirna_data_clean != 0).any(axis=1)]
    
    print(f"Final dataset: {mirna_data_clean.shape[0]} miRNAs, {mirna_data_clean.shape[1]} samples")
    
    # Save cleaned data
    output_file = '../data/mirna_counts_cleaned.csv'
    mirna_data_clean.to_csv(output_file)
    print(f"Cleaned miRNA data saved to {output_file}")
    
    return mirna_data_clean

def preprocess_methylation_data():
    """Preprocess DNA methylation data."""
    print("Preprocessing DNA methylation data...")
    
    methylation_data = pd.read_csv('../data/merged-WGBS-CpG-counts_filtered.csv')
    
    # Clean column names - extract just the ACR-XXX-TPX part
    clean_cols = {}
    for col in methylation_data.columns:
        if 'ACR-' in col and 'TP' in col:
            # Extract ACR-XXX-TPX pattern
            match = re.search(r'ACR-\d+-TP\d+', col)
            if match:
                clean_name = match.group()
                clean_cols[col] = clean_name
                print(f"  {col} -> {clean_name}")
    
    # Rename columns
    methylation_data = methylation_data.rename(columns=clean_cols)
    
    # Keep only essential columns
    essential_cols = ['CpG'] + list(clean_cols.values())
    methylation_data_clean = methylation_data[essential_cols].copy()
    
    # Set CpG as index
    methylation_data_clean.set_index('CpG', inplace=True)
    
    # Remove rows with all zero values
    methylation_data_clean = methylation_data_clean.loc[(methylation_data_clean != 0).any(axis=1)]
    
    # Fill NaN values with 0
    methylation_data_clean = methylation_data_clean.fillna(0)
    
    print(f"Final dataset: {methylation_data_clean.shape[0]} CpGs, {methylation_data_clean.shape[1]} samples")
    
    # Save cleaned data
    output_file = '../data/methylation_counts_cleaned.csv'
    methylation_data_clean.to_csv(output_file)
    print(f"Cleaned methylation data saved to {output_file}")
    
    return methylation_data_clean

def preprocess_gene_data():
    """Preprocess gene expression data."""
    print("Preprocessing gene expression data...")
    
    gene_data = pd.read_csv('../data/apul-gene_count_matrix.csv')
    
    # Remove duplicate columns (some samples appear twice)
    gene_data = gene_data.loc[:, ~gene_data.columns.duplicated()]
    
    # Clean column names - extract just the ACR-XXX-TPX part
    clean_cols = {}
    for col in gene_data.columns:
        if col == 'gene_id':
            continue
        if 'ACR-' in col and 'TP' in col:
            # Extract ACR-XXX-TPX pattern
            match = re.search(r'ACR-\d+-TP\d+', col)
            if match:
                clean_name = match.group()
                clean_cols[col] = clean_name
                print(f"  {col} -> {clean_name}")
    
    # Rename columns
    gene_data = gene_data.rename(columns=clean_cols)
    
    # Keep only essential columns
    essential_cols = ['gene_id'] + list(clean_cols.values())
    gene_data_clean = gene_data[essential_cols].copy()
    
    # Set gene_id as index
    gene_data_clean.set_index('gene_id', inplace=True)
    
    # Remove genes with all zero counts
    gene_data_clean = gene_data_clean.loc[(gene_data_clean != 0).any(axis=1)]
    
    print(f"Final dataset: {gene_data_clean.shape[0]} genes, {gene_data_clean.shape[1]} samples")
    
    # Save cleaned data
    output_file = '../data/gene_counts_cleaned.csv'
    gene_data_clean.to_csv(output_file)
    print(f"Cleaned gene data saved to {output_file}")
    
    return gene_data_clean

def main():
    """Main preprocessing pipeline."""
    print("Data Preprocessing Pipeline")
    print("=" * 40)
    
    # Preprocess all data
    gene_data = preprocess_gene_data()
    print()
    
    mirna_data = preprocess_mirna_data()
    print()
    
    methylation_data = preprocess_methylation_data()
    print()
    
    # Check sample overlap
    print("Sample Overlap Analysis:")
    print("-" * 30)
    
    gene_samples = set(gene_data.columns)
    mirna_samples = set(mirna_data.columns)
    methylation_samples = set(methylation_data.columns)
    
    print(f"Gene samples: {len(gene_samples)}")
    print(f"miRNA samples: {len(mirna_samples)}")
    print(f"Methylation samples: {len(methylation_samples)}")
    
    # Find common samples
    common_samples = gene_samples & mirna_samples & methylation_samples
    print(f"Common samples across all datasets: {len(common_samples)}")
    
    if len(common_samples) > 0:
        print("✓ Good sample overlap detected")
        print(f"Common samples: {sorted(list(common_samples))[:10]}...")
    else:
        print("⚠ Warning: No common samples found across all datasets")
    
    print("\nPreprocessing complete!")

if __name__ == "__main__":
    main()
