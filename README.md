# Timeseries-02: Multi-Omics Gene Regulation Analysis Pipeline

A comprehensive bioinformatics pipeline for analyzing the influence of miRNA, lncRNA, and DNA methylation on gene expression using multi-omics data from 40 samples across multiple time points and conditions.

## 🧬 Project Overview

This repository contains a sophisticated analysis pipeline that examines how three regulatory layers influence gene expression patterns:
- **miRNA Expression** - microRNA regulatory molecules
- **lncRNA Expression** - long non-coding RNA molecules  
- **DNA Methylation** - epigenetic modifications at CpG sites

The pipeline applies multiple machine learning approaches to quantify regulatory influence and identify key regulatory relationships across time series data.

## 📁 Repository Structure

```
timeseries-02/
├── code/                          # Analysis scripts and pipeline
│   ├── comprehensive_regulation_analysis.py  # Main integrated analysis script
│   ├── clean_and_standardize_datasets.py     # Data cleaning and standardization
│   ├── verify_cleaned_datasets.py            # Data validation scripts
│   ├── analyze_gene_regulation.py            # Gene regulation analysis
│   ├── simple_correlation_analysis.py        # Correlation-based analysis
│   ├── preprocess_all_data.py                # Data preprocessing pipeline
│   ├── preprocess_lncrna.py                  # lncRNA-specific preprocessing
│   ├── test_data_loading.py                  # Data loading tests
│   ├── run_analysis.sh                       # Automated pipeline runner
│   ├── requirements.txt                      # Python dependencies
│   └── README.md                             # Detailed code documentation
├── data/                          # Input data files
│   ├── cleaned_datasets/                     # Processed and cleaned data
│   │   ├── gene_counts_cleaned.csv          # Cleaned gene expression data
│   │   ├── lncrna_counts_cleaned.csv        # Cleaned lncRNA expression data
│   │   ├── mirna_counts_cleaned.csv         # Cleaned miRNA expression data
│   │   └── wgbs_counts_cleaned.csv          # Cleaned DNA methylation data
│   ├── apul-gene_count_matrix.csv           # Raw gene expression (RNA-seq)
│   ├── Apul_miRNA_counts_formatted.txt     # Raw miRNA expression
│   ├── Apul_lncRNA_counts_filtered.txt     # Raw lncRNA expression
│   └── merged-WGBS-CpG-counts_filtered.csv # Raw DNA methylation
├── output/                        # Analysis results and visualizations
│   ├── integrated_dataset.csv               # Combined multi-omics dataset
│   ├── regulatory_influence_results.csv     # Model performance results
│   ├── correlation_results.csv              # Correlation analysis results
│   ├── analysis_summary.txt                 # Analysis summary report
│   ├── correlation_summary.txt              # Correlation summary
│   ├── model_comparison.png                 # Model performance comparison
│   ├── r2_distribution.png                  # R² score distributions
│   ├── sample_correlation.png               # Sample correlation heatmap
│   ├── significant_correlations.png         # Significant correlations
│   ├── correlation_distributions.png        # Correlation distributions
│   └── correlation_heatmap.png             # Correlation heatmap
└── README.md                      # This file
```

## 🔬 Analysis Methods

### 1. Data Preprocessing Pipeline
- **Data Cleaning**: Multi-format data loading and cleaning
- **Sample Standardization**: Sample name standardization across datasets
- **Quality Filtering**: Zero-expression removal and data validation
- **Data Integration**: Harmonized multi-omics dataset creation

### 2. Comprehensive Analysis Approaches
- **Correlation Analysis**: Cross-dataset correlations between regulatory layers
- **Time Series Analysis**: Temporal patterns across multiple time points
- **Regulatory Network Inference**: Identification of regulatory relationships
- **Machine Learning Models**: Multiple approaches for regulatory influence quantification

### 3. Machine Learning Models
- **Multiple Linear Regression**: Basic linear relationships
- **Ridge Regression**: Regularized linear regression
- **Lasso Regression**: Sparse regression with feature selection
- **Random Forest**: Non-linear ensemble method
- **Cross-Validation**: Robust model evaluation with 3-fold CV

### 4. Statistical Evaluation
- **R² Scores**: Model performance metrics for each gene
- **Correlation Coefficients**: Pearson and Spearman correlations
- **Significance Testing**: Statistical validation of relationships
- **Feature Importance**: Regulatory factor contribution analysis

## 🚀 Quick Start

### Prerequisites
- Python 3.7+
- Required packages (see `code/requirements.txt`)

### Installation
```bash
# Clone the repository
git clone <repository-url>
cd timeseries-02

# Install dependencies
pip install -r code/requirements.txt
```

### Run Analysis
```bash
# Automated pipeline (recommended)
bash code/run_analysis.sh

# Or run individual components
cd code
python comprehensive_regulation_analysis.py    # Main analysis
python simple_correlation_analysis.py         # Correlation analysis
python clean_and_standardize_datasets.py     # Data preprocessing
```

## 📊 Output Files

The pipeline generates comprehensive results including:

### Data Files
- **`integrated_dataset.csv`**: Combined multi-omics dataset ready for analysis
- **`regulatory_influence_results.csv`**: Model performance results for each gene
- **`correlation_results.csv`**: Correlation analysis results between regulatory layers

### Summary Reports
- **`analysis_summary.txt`**: Overall analysis summary and key findings
- **`correlation_summary.txt`**: Summary of correlation analysis results

### Visualizations
- **`model_comparison.png`**: Performance comparison across different ML models
- **`r2_distribution.png`**: Distribution of R² scores across genes
- **`sample_correlation.png`**: Sample-to-sample correlation heatmap
- **`correlation_heatmap.png`**: Regulatory factor correlation matrix
- **`significant_correlations.png`**: Visualization of significant correlations
- **`correlation_distributions.png`**: Distribution of correlation coefficients

## 🧪 Data Sources

- **Gene Expression**: RNA-seq count matrix from 40 samples across multiple time points
- **miRNA**: Small RNA sequencing data (cleaned and formatted)
- **lncRNA**: Long non-coding RNA expression profiles (filtered and processed)
- **DNA Methylation**: Whole genome bisulfite sequencing (WGBS) data

## 🔍 Interpretation Guide

### R² Score Meaning
- **R² = 1.0**: Perfect prediction (100% variance explained)
- **R² = 0.5**: 50% variance explained by regulatory factors
- **R² = 0.0**: No predictive power
- **R² < 0**: Model performs worse than baseline

### Regulatory Influence Levels
- **High (R² > 0.5)**: Strong regulatory influence
- **Medium (0.2 < R² < 0.5)**: Moderate regulatory influence
- **Low (R² < 0.2)**: Weak regulatory influence

### Correlation Interpretation
- **Strong Positive (> 0.7)**: Strong positive regulatory relationship
- **Moderate Positive (0.3-0.7)**: Moderate positive relationship
- **Weak (-0.3 to 0.3)**: Weak or no relationship
- **Moderate Negative (-0.7 to -0.3)**: Moderate negative relationship
- **Strong Negative (< -0.7)**: Strong negative regulatory relationship

## 🛠️ Technical Details

- **Cross-validation**: 3-fold CV for robust evaluation
- **Feature scaling**: StandardScaler for regulatory factors
- **Data integration**: Averaging across molecular features per sample
- **Sample matching**: Automatic sample name standardization
- **Memory optimization**: Efficient handling of large datasets
- **Time series support**: Analysis across multiple time points and conditions

## 📚 Dependencies

- **pandas** ≥1.3.0: Data manipulation and analysis
- **numpy** ≥1.21.0: Numerical computing
- **matplotlib** ≥3.4.0: Plotting and visualization
- **seaborn** ≥0.11.0: Statistical data visualization
- **scipy** ≥1.7.0: Scientific computing
- **scikit-learn** ≥1.0.0: Machine learning algorithms

## 🔄 Workflow

1. **Data Preprocessing**: Clean and standardize raw multi-omics data
2. **Data Integration**: Create harmonized dataset with aligned samples
3. **Correlation Analysis**: Examine relationships between regulatory layers
4. **Machine Learning**: Apply multiple models to quantify regulatory influence
5. **Statistical Evaluation**: Assess model performance and significance
6. **Visualization**: Generate comprehensive plots and summaries
7. **Results Export**: Save analysis results and visualizations

## 🤝 Contributing

This pipeline is designed for bioinformatics research and can be extended with:
- Additional regulatory layers (e.g., transcription factors, chromatin accessibility)
- New machine learning models (e.g., deep learning, graph neural networks)
- Enhanced visualization options and interactive plots
- Integration with other omics data types (e.g., proteomics, metabolomics)
- Time series analysis improvements and temporal modeling

## 📄 License

[Add your license information here]

## 👥 Authors

[Add author information here]

## 📖 Citation

If you use this pipeline in your research, please cite:
[Add citation information here]

## 📞 Support

For questions or issues with the pipeline, please:
1. Check the `code/README.md` for detailed technical documentation
2. Review the output files for analysis results
3. Examine the code comments for implementation details