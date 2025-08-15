#!/usr/bin/env python3
"""
Script to clean and standardize all datasets for correlation analysis.
This script:
1. Creates a subdirectory for cleaned datasets
2. Standardizes sample IDs across all datasets
3. Removes genes with zero expression across all samples
4. Removes genes with limited variation across samples
5. Ensures common structure and sample IDs
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def setup_directories():
    """Create necessary directories for cleaned datasets."""
    cleaned_dir = Path("data/cleaned_datasets")
    cleaned_dir.mkdir(exist_ok=True)
    logger.info(f"Created directory: {cleaned_dir}")
    return cleaned_dir

def load_gene_counts():
    """Load and clean gene count matrix."""
    logger.info("Loading gene count matrix...")
    
    # Load the large gene count file
    gene_counts = pd.read_csv("data/apul-gene_count_matrix.csv", index_col=0)
    
    # Clean column names - remove any whitespace and duplicates
    gene_counts.columns = gene_counts.columns.str.strip()
    
    # Remove duplicate columns if they exist
    gene_counts = gene_counts.loc[:, ~gene_counts.columns.duplicated()]
    
    logger.info(f"Loaded {gene_counts.shape[0]} genes and {gene_counts.shape[1]} samples")
    return gene_counts

def load_lncrna_counts():
    """Load and clean lncRNA count matrix."""
    logger.info("Loading lncRNA count matrix...")
    
    # Use the cleaned version if available
    lncrna_file = "data/lncrna_counts_cleaned.csv"
    if os.path.exists(lncrna_file):
        lncrna_counts = pd.read_csv(lncrna_file, index_col=0)
    else:
        # Fall back to the filtered version
        lncrna_counts = pd.read_csv("data/Apul_lncRNA_counts_filtered.txt", 
                                   sep='\t', comment='#', index_col=0)
        # Extract only the count columns (skip metadata columns)
        count_cols = [col for col in lncrna_counts.columns if 'ACR-' in col]
        lncrna_counts = lncrna_counts[count_cols]
    
    logger.info(f"Loaded {lncrna_counts.shape[0]} lncRNAs and {lncrna_counts.shape[1]} samples")
    return lncrna_counts

def load_mirna_counts():
    """Load and clean miRNA count matrix."""
    logger.info("Loading miRNA count matrix...")
    
    # Load miRNA counts
    mirna_counts = pd.read_csv("data/Apul_miRNA_counts_formatted.txt", sep='\t', index_col=0)
    
    # Clean column names - convert to standard ACR format
    # Current format: 1A10_ACR-145_TP4 -> ACR-145-TP4
    new_cols = []
    for col in mirna_counts.columns:
        if 'ACR-' in col:
            # Extract ACR and TP information
            parts = col.split('_')
            for part in parts:
                if 'ACR-' in part:
                    acr_part = part
                elif 'TP' in part:
                    tp_part = part
            new_col = f"{acr_part}-{tp_part}"
            new_cols.append(new_col)
        else:
            new_cols.append(col)
    
    mirna_counts.columns = new_cols
    
    logger.info(f"Loaded {mirna_counts.shape[0]} miRNAs and {mirna_counts.shape[1]} samples")
    return mirna_counts

def load_wgbs_counts():
    """Load and clean WGBS CpG count matrix."""
    logger.info("Loading WGBS CpG count matrix...")
    
    wgbs_counts = pd.read_csv("data/merged-WGBS-CpG-counts_filtered.csv", index_col=0)
    
    # Clean column names - extract ACR-TP format
    new_cols = []
    for col in wgbs_counts.columns:
        if 'ACR-' in col:
            # Extract ACR-TP format from the long filename
            parts = col.split('.')
            for part in parts:
                if 'ACR-' in part and 'TP' in part:
                    new_cols.append(part)
                    break
            else:
                new_cols.append(col)
        else:
            new_cols.append(col)
    
    wgbs_counts.columns = new_cols
    
    logger.info(f"Loaded {wgbs_counts.shape[0]} CpGs and {wgbs_counts.shape[1]} samples")
    return wgbs_counts

def standardize_sample_ids(df, dataset_name):
    """Standardize sample IDs to ACR-XXX-TPX format."""
    logger.info(f"Standardizing sample IDs for {dataset_name}...")
    
    # Define the expected sample format: ACR-XXX-TPX
    expected_samples = []
    for acr in [139, 145, 150, 173, 186, 225, 229, 237, 244, 265]:
        for tp in [1, 2, 3, 4]:
            expected_samples.append(f"ACR-{acr}-TP{tp}")
    
    # Find which expected samples are present in the dataset
    available_samples = []
    for expected in expected_samples:
        # Try exact match first
        if expected in df.columns:
            available_samples.append(expected)
        else:
            # Try partial match
            matching_cols = [col for col in df.columns if expected in col]
            if matching_cols:
                available_samples.append(matching_cols[0])
    
    logger.info(f"Found {len(available_samples)} matching samples out of {len(expected_samples)} expected")
    
    # Keep only the available samples
    df_standardized = df[available_samples].copy()
    
    # Rename columns to standard format if needed
    for i, col in enumerate(available_samples):
        if col != expected_samples[i]:
            df_standardized = df_standardized.rename(columns={col: expected_samples[i]})
    
    return df_standardized

def remove_zero_expression(df, dataset_name):
    """Remove features with zero expression across all samples."""
    logger.info(f"Removing zero expression features from {dataset_name}...")
    
    # Count non-zero values per row
    non_zero_counts = (df != 0).sum(axis=1)
    
    # Keep only features with at least one non-zero value
    df_filtered = df[non_zero_counts > 0].copy()
    
    removed_count = len(df) - len(df_filtered)
    logger.info(f"Removed {removed_count} features with zero expression across all samples")
    
    return df_filtered

def remove_low_variation(df, dataset_name, min_cv=0.1):
    """Remove features with limited variation across samples (low coefficient of variation)."""
    logger.info(f"Removing low variation features from {dataset_name} (min CV: {min_cv})...")
    
    # Calculate coefficient of variation for each feature
    # CV = std / mean, but handle cases where mean is 0
    means = df.mean(axis=1)
    stds = df.std(axis=1)
    
    # Calculate CV, setting CV to 0 where mean is 0
    cv = np.where(means > 0, stds / means, 0)
    
    # Keep features with CV above threshold
    high_var_mask = cv >= min_cv
    df_filtered = df[high_var_mask].copy()
    
    removed_count = len(df) - len(df_filtered)
    logger.info(f"Removed {removed_count} features with low variation (CV < {min_cv})")
    
    return df_filtered

def ensure_common_samples(datasets):
    """Ensure all datasets have the same sample IDs."""
    logger.info("Ensuring common sample IDs across all datasets...")
    
    # Get the intersection of all sample IDs
    dataset_names = list(datasets.keys())
    common_samples = set(datasets[dataset_names[0]].columns)
    for name, df in datasets.items():
        common_samples = common_samples.intersection(set(df.columns))
    
    common_samples = sorted(list(common_samples))
    logger.info(f"Found {len(common_samples)} common samples across all datasets")
    
    # Filter all datasets to common samples
    datasets_common = {}
    for name, df in datasets.items():
        datasets_common[name] = df[common_samples].copy()
        logger.info(f"{name}: {datasets_common[name].shape}")
    
    return datasets_common

def save_cleaned_datasets(datasets, output_dir):
    """Save all cleaned datasets."""
    logger.info("Saving cleaned datasets...")
    
    for name, df in datasets.items():
        output_file = output_dir / f"{name}_cleaned.csv"
        df.to_csv(output_file)
        logger.info(f"Saved {name} to {output_file}")
        
        # Also save a summary
        summary_file = output_dir / f"{name}_summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"Dataset: {name}\n")
            f.write(f"Features: {df.shape[0]}\n")
            f.write(f"Samples: {df.shape[1]}\n")
            f.write(f"Sample IDs: {', '.join(df.columns)}\n")
            f.write(f"Total non-zero values: {(df != 0).sum().sum()}\n")
            f.write(f"Mean expression per feature: {df.mean(axis=1).mean():.2f}\n")
            f.write(f"Std expression per feature: {df.std(axis=1).mean():.2f}\n")

def main():
    """Main function to clean and standardize all datasets."""
    logger.info("Starting dataset cleaning and standardization...")
    
    # Setup directories
    output_dir = setup_directories()
    
    try:
        # Load all datasets
        gene_counts = load_gene_counts()
        lncrna_counts = load_lncrna_counts()
        mirna_counts = load_mirna_counts()
        wgbs_counts = load_wgbs_counts()
        
        # Standardize sample IDs
        gene_counts_std = standardize_sample_ids(gene_counts, "Gene counts")
        lncrna_counts_std = standardize_sample_ids(lncrna_counts, "lncRNA counts")
        mirna_counts_std = standardize_sample_ids(mirna_counts, "miRNA counts")
        wgbs_counts_std = standardize_sample_ids(wgbs_counts, "WGBS CpG counts")
        
        # Remove zero expression features
        gene_counts_clean = remove_zero_expression(gene_counts_std, "Gene counts")
        lncrna_counts_clean = remove_zero_expression(lncrna_counts_std, "lncRNA counts")
        mirna_counts_clean = remove_zero_expression(mirna_counts_std, "miRNA counts")
        wgbs_counts_clean = remove_zero_expression(wgbs_counts_std, "WGBS CpG counts")
        
        # Remove low variation features
        gene_counts_final = remove_low_variation(gene_counts_clean, "Gene counts", min_cv=0.1)
        lncrna_counts_final = remove_low_variation(lncrna_counts_clean, "lncRNA counts", min_cv=0.1)
        mirna_counts_final = remove_low_variation(mirna_counts_clean, "miRNA counts", min_cv=0.1)
        wgbs_counts_final = remove_low_variation(wgbs_counts_clean, "WGBS CpG counts", min_cv=0.1)
        
        # Ensure common samples across all datasets
        datasets = {
            "gene_counts": gene_counts_final,
            "lncrna_counts": lncrna_counts_final,
            "mirna_counts": mirna_counts_final,
            "wgbs_counts": wgbs_counts_final
        }
        
        datasets_common = ensure_common_samples(datasets)
        
        # Save cleaned datasets
        save_cleaned_datasets(datasets_common, output_dir)
        
        # Create a combined summary
        summary_file = output_dir / "combined_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("COMBINED DATASET SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Common samples: {len(datasets_common[list(datasets_common.keys())[0]].columns)}\n")
            f.write(f"Sample IDs: {', '.join(datasets_common[list(datasets_common.keys())[0]].columns)}\n\n")
            
            for name, df in datasets_common.items():
                f.write(f"{name.upper()}:\n")
                f.write(f"  Features: {df.shape[0]}\n")
                f.write(f"  Samples: {df.shape[1]}\n")
                f.write(f"  Sparsity: {((df == 0).sum().sum() / (df.shape[0] * df.shape[1]) * 100):.1f}%\n\n")
        
        logger.info("Dataset cleaning and standardization completed successfully!")
        logger.info(f"All cleaned datasets saved to: {output_dir}")
        
    except Exception as e:
        logger.error(f"Error during dataset cleaning: {str(e)}")
        raise

if __name__ == "__main__":
    main()
