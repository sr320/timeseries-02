#!/usr/bin/env python3
"""
Verification script to confirm that all cleaned datasets have the same structure
and sample IDs for correlation analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path

def verify_datasets():
    """Verify that all cleaned datasets have consistent structure."""
    
    # Load all cleaned datasets
    cleaned_dir = Path("data/cleaned_datasets")
    
    datasets = {}
    for file in cleaned_dir.glob("*_cleaned.csv"):
        name = file.stem.replace("_cleaned", "")
        datasets[name] = pd.read_csv(file, index_col=0)
        print(f"Loaded {name}: {datasets[name].shape}")
    
    print("\n" + "="*60)
    print("DATASET VERIFICATION RESULTS")
    print("="*60)
    
    # Check sample consistency
    sample_sets = [set(df.columns) for df in datasets.values()]
    all_samples_match = all(sample_sets[0] == s for s in sample_sets)
    
    if all_samples_match:
        print("✓ All datasets have identical sample IDs")
        print(f"  Common samples: {len(sample_sets[0])}")
        print(f"  Sample IDs: {', '.join(sorted(sample_sets[0]))}")
    else:
        print("✗ Sample ID mismatch detected!")
        for i, (name, samples) in enumerate(zip(datasets.keys(), sample_sets)):
            print(f"  {name}: {len(samples)} samples")
    
    print("\n" + "-"*60)
    
    # Check for zero expression across all samples
    print("ZERO EXPRESSION CHECK:")
    for name, df in datasets.items():
        zero_rows = (df == 0).all(axis=1).sum()
        total_rows = len(df)
        print(f"  {name}: {zero_rows}/{total_rows} features with zero expression across all samples")
    
    print("\n" + "-"*60)
    
    # Check variation statistics
    print("VARIATION STATISTICS:")
    for name, df in datasets.items():
        # Calculate coefficient of variation for each feature
        means = df.mean(axis=1)
        stds = df.std(axis=1)
        cv = np.where(means > 0, stds / means, 0)
        
        min_cv = cv[cv > 0].min() if (cv > 0).any() else 0
        max_cv = cv.max()
        mean_cv = cv[cv > 0].mean() if (cv > 0).any() else 0
        
        print(f"  {name}:")
        print(f"    Min CV: {min_cv:.3f}")
        print(f"    Max CV: {max_cv:.3f}")
        print(f"    Mean CV: {mean_cv:.3f}")
    
    print("\n" + "-"*60)
    
    # Check data types and missing values
    print("DATA QUALITY CHECK:")
    for name, df in datasets.items():
        missing_values = df.isnull().sum().sum()
        data_types = df.dtypes.value_counts()
        print(f"  {name}:")
        print(f"    Missing values: {missing_values}")
        print(f"    Data types: {dict(data_types)}")
    
    print("\n" + "="*60)
    
    if all_samples_match:
        print("✓ VERIFICATION PASSED: All datasets are ready for correlation analysis!")
    else:
        print("✗ VERIFICATION FAILED: Sample ID inconsistencies detected!")
    
    return all_samples_match

if __name__ == "__main__":
    verify_datasets()
