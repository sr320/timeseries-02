# Gene Regulation Analysis Pipeline

This pipeline analyzes the influence of miRNA, lncRNA, and DNA methylation on gene expression using multi-omics data from 40 samples.

## Overview

The analysis examines how three regulatory layers (miRNA, lncRNA, and DNA methylation) influence gene expression patterns. Using a subset of 20 highly variable genes, the pipeline applies multiple regression models to quantify the regulatory influence of each molecular layer.

## Data Input

- **Gene Expression**: `data/apul-gene_count_matrix.csv` - RNA-seq count matrix
- **miRNA Expression**: `data/Apul_miRNA_counts_formatted.txt` - miRNA count matrix  
- **lncRNA Expression**: `data/Apul_lncRNA_counts_filtered.txt` - lncRNA count matrix
- **DNA Methylation**: `data/merged-WGBS-CpG-counts_filtered.csv` - CpG methylation levels

## Analysis Pipeline

### 1. Data Preprocessing
- Load and clean all data matrices
- Standardize sample names across datasets
- Apply log2 transformation to count data
- Remove genes/features with zero expression

### 2. Gene Selection
- Select top 20 genes by variance (most variable genes)
- These genes are likely to show strong regulatory patterns

### 3. Regulatory Influence Analysis
- **Multiple Linear Regression**: Basic linear relationship
- **Ridge Regression**: Regularized linear regression
- **Lasso Regression**: Sparse linear regression with feature selection
- **Random Forest**: Non-linear ensemble method

### 4. Output Generation
- Statistical results for each gene and model
- Comprehensive visualizations
- Summary statistics and rankings

## Usage

### Install Dependencies
```bash
pip install -r requirements.txt
```

### Run Analysis
```bash
python analyze_gene_regulation.py
```

## Output Files

The analysis generates several output files in the `output/` directory:

- **`regulatory_influence_results.csv`**: Detailed R² scores for each gene and model
- **`integrated_dataset.csv`**: Combined dataset with all molecular layers
- **`analysis_summary.txt`**: Summary statistics and top-performing genes
- **`regulatory_influence_heatmap.png`**: Heatmap showing model performance by gene
- **`r2_distribution.png`**: Distribution of regulatory influence scores
- **`model_comparison.png`**: Boxplot comparing model performance
- **`sample_correlation.png`**: Sample correlation matrix

## Interpretation

### R² Score Meaning
- **R² = 1.0**: Perfect prediction (100% of variance explained)
- **R² = 0.5**: 50% of variance explained by regulatory factors
- **R² = 0.0**: No predictive power
- **R² < 0**: Model performs worse than baseline

### Regulatory Influence Levels
- **High (R² > 0.5)**: Strong regulatory influence
- **Medium (0.2 < R² < 0.5)**: Moderate regulatory influence  
- **Low (R² < 0.2)**: Weak regulatory influence

## Technical Details

- **Cross-validation**: 3-fold cross-validation for robust model evaluation
- **Feature scaling**: StandardScaler applied to regulatory factors
- **Data integration**: Averaging across miRNAs, lncRNAs, and CpGs per sample
- **Sample matching**: Automatic sample name standardization across datasets

## Troubleshooting

- Ensure all input files are in the `data/` directory
- Check that sample names contain "ACR-" and "TP" patterns
- Verify sufficient memory for large datasets
- Check Python package versions match requirements.txt
