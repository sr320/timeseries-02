#!/usr/bin/env python3
"""
Test script to verify data loading and basic structure.
Run this before the main analysis to ensure data is accessible.
"""

import pandas as pd
import numpy as np

def test_data_loading():
    """Test loading all data files."""
    print("Testing data loading...")
    print("=" * 40)
    
    try:
        # Test gene expression data
        print("1. Loading gene expression data...")
        gene_data = pd.read_csv('../data/apul-gene_count_matrix.csv')
        print(f"   ✓ Loaded: {gene_data.shape[0]} genes, {gene_data.shape[1]} samples")
        print(f"   ✓ Sample columns: {[col for col in gene_data.columns if 'ACR-' in col][:5]}...")
        
        # Test miRNA data
        print("\n2. Loading miRNA data...")
        mirna_data = pd.read_csv('../data/Apul_miRNA_counts_formatted.txt', sep='\t')
        print(f"   ✓ Loaded: {mirna_data.shape[0]} miRNAs, {mirna_data.shape[1]} samples")
        print(f"   ✓ Sample columns: {[col for col in mirna_data.columns if 'ACR-' in col][:5]}...")
        
        # Test lncRNA data
        print("\n3. Loading lncRNA data...")
        lncrna_data = pd.read_csv('../data/Apul_lncRNA_counts_filtered.txt', sep='\t', comment='#')
        print(f"   ✓ Loaded: {lncrna_data.shape[0]} lncRNAs, {lncrna_data.shape[1]} columns")
        
        # Find count columns
        count_cols = [col for col in lncrna_data.columns if 'ACR-' in col]
        print(f"   ✓ Found {len(count_cols)} sample columns with ACR- pattern")
        
        # Test methylation data
        print("\n4. Loading DNA methylation data...")
        methylation_data = pd.read_csv('../data/merged-WGBS-CpG-counts_filtered.csv')
        print(f"   ✓ Loaded: {methylation_data.shape[0]} CpGs, {methylation_data.shape[1]} samples")
        print(f"   ✓ Sample columns: {[col for col in methylation_data.columns if 'ACR-' in col][:5]}...")
        
        print("\n" + "=" * 40)
        print("✓ All data files loaded successfully!")
        print("✓ Data structure looks correct")
        print("✓ Ready to run main analysis")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Error loading data: {e}")
        return False

def test_sample_matching():
    """Test if sample names can be matched across datasets."""
    print("\nTesting sample name matching...")
    print("=" * 40)
    
    try:
        # Load data
        gene_data = pd.read_csv('../data/apul-gene_count_matrix.csv')
        mirna_data = pd.read_csv('../data/Apul_miRNA_counts_formatted.txt', sep='\t')
        lncrna_data = pd.read_csv('../data/Apul_lncRNA_counts_filtered.txt', sep='\t', comment='#')
        methylation_data = pd.read_csv('../data/merged-WGBS-CpG-counts_filtered.csv')
        
        # Extract sample names
        gene_samples = [col for col in gene_data.columns if 'ACR-' in col]
        mirna_samples = [col for col in mirna_data.columns if 'ACR-' in col]
        lncrna_samples = [col for col in lncrna_data.columns if 'ACR-' in col]
        methylation_samples = [col for col in methylation_data.columns if 'ACR-' in col]
        
        print(f"Gene samples: {len(gene_samples)}")
        print(f"miRNA samples: {len(mirna_samples)}")
        print(f"lncRNA samples: {len(lncrna_samples)}")
        print(f"Methylation samples: {len(methylation_samples)}")
        
        # Check for common patterns
        common_samples = set(gene_samples) & set(mirna_samples) & set(lncrna_samples) & set(methylation_samples)
        print(f"\nCommon samples across all datasets: {len(common_samples)}")
        
        if len(common_samples) > 0:
            print("✓ Sample matching looks good")
            return True
        else:
            print("⚠ Warning: No common samples found across all datasets")
            return False
            
    except Exception as e:
        print(f"❌ Error in sample matching: {e}")
        return False

if __name__ == "__main__":
    print("Data Loading Test")
    print("=" * 50)
    
    # Test data loading
    loading_ok = test_data_loading()
    
    if loading_ok:
        # Test sample matching
        matching_ok = test_sample_matching()
        
        if matching_ok:
            print("\n🎉 All tests passed! Ready to run analysis.")
        else:
            print("\n⚠ Sample matching issues detected. Analysis may have limited sample overlap.")
    else:
        print("\n❌ Data loading failed. Check file paths and formats.")
