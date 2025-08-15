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
│   │   ├── wgbs_counts_cleaned.csv          # Cleaned DNA methylation data
│   │   ├── gene_counts_summary.txt          # Gene data summary statistics
│   │   ├── lncrna_counts_summary.txt        # lncRNA data summary statistics
│   │   ├── mirna_counts_summary.txt         # miRNA data summary statistics
│   │   ├── wgbs_counts_summary.txt          # Methylation data summary statistics
│   │   ├── combined_summary.txt             # Overall dataset summary
│   │   └── README.md                        # Data documentation
│   ├── apul-gene_count_matrix.csv           # Raw gene expression (RNA-seq)
│   ├── Apul_miRNA_counts_formatted.txt     # Raw miRNA expression
│   ├── Apul_lncRNA_counts_filtered.txt     # Raw lncRNA expression
│   └── merged-WGBS-CpG-counts_filtered.csv # Raw DNA methylation
├── output/                        # Analysis results and visualizations
│   ├── regulation_analysis/                # Comprehensive regulation analysis results
│   │   ├── COMPREHENSIVE_ANALYSIS_SUMMARY_REPORT.md  # Detailed analysis report
│   │   ├── model_performance_summary.csv             # Model performance summary
│   │   ├── time_series_results.json                  # Time series analysis results
│   │   ├── methylation_gene_network.csv              # Methylation-gene regulatory network
│   │   ├── lncrna_gene_network.csv                  # lncRNA-gene regulatory network
│   │   ├── mirna_gene_network.csv                   # miRNA-gene regulatory network
│   │   ├── lncrna_mirna_feature_correlations.csv    # lncRNA-miRNA feature correlations
│   │   ├── gene_methylation_feature_correlations.csv # Gene-methylation feature correlations
│   │   ├── gene_mirna_feature_correlations.csv      # Gene-miRNA feature correlations
│   │   ├── gene_lncrna_feature_correlations.csv     # Gene-lncRNA feature correlations
│   │   ├── lncrna_mirna_sample_correlations.csv     # lncRNA-miRNA sample correlations
│   │   ├── gene_methylation_sample_correlations.csv  # Gene-methylation sample correlations
│   │   ├── gene_mirna_sample_correlations.csv       # Gene-miRNA sample correlations
│   │   ├── gene_lncrna_sample_correlations.csv      # Gene-lncRNA sample correlations
│   │   ├── model_performance_comparison.png         # Model performance comparison
│   │   ├── regulation_type_distributions.png        # Regulation type distributions
│   │   ├── time_series_analysis.png                 # Time series analysis plots
│   │   └── sample_correlation_heatmaps.png         # Sample correlation heatmaps
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
- **Summary Statistics**: Comprehensive data quality reports for each dataset

### 2. Comprehensive Analysis Approaches
- **Correlation Analysis**: Cross-dataset correlations between regulatory layers
- **Time Series Analysis**: Temporal patterns across multiple time points (TP1-TP4)
- **Regulatory Network Inference**: Identification of regulatory relationships
- **Machine Learning Models**: Multiple approaches for regulatory influence quantification
- **Network Analysis**: Detailed regulatory network construction and analysis

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
- **Network Metrics**: Regulatory relationship density and patterns

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

### Comprehensive Regulation Analysis
- **`COMPREHENSIVE_ANALYSIS_SUMMARY_REPORT.md`**: Detailed analysis report with biological insights
- **`model_performance_summary.csv`**: Comprehensive model performance metrics
- **`time_series_results.json`**: Time series analysis results across conditions
- **Regulatory Network Files**: Detailed regulatory relationship matrices for each data type combination
- **Feature Correlation Files**: Comprehensive correlation analysis between regulatory layers

### Visualizations
- **`model_comparison.png`**: Performance comparison across different ML models
- **`r2_distribution.png`**: Distribution of R² scores across genes
- **`sample_correlation.png`**: Sample-to-sample correlation heatmap
- **`correlation_heatmap.png`**: Regulatory factor correlation matrix
- **`significant_correlations.png`**: Visualization of significant correlations
- **`correlation_distributions.png`**: Distribution of correlation coefficients
- **`regulation_type_distributions.png`**: Distribution of regulation types
- **`time_series_analysis.png`**: Temporal dynamics across conditions
- **`sample_correlation_heatmaps.png`**: Comprehensive sample correlation analysis

## 🧪 Data Sources

- **Gene Expression**: RNA-seq count matrix from 40 samples across 10 conditions × 4 time points
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

### Regulatory Network Insights
- **miRNA-Gene Networks**: 3,446 regulatory relationships identified
- **lncRNA-Gene Networks**: 40,652 regulatory relationships identified
- **Methylation-Gene Networks**: 9,214 regulatory relationships identified

## 🛠️ Technical Details

- **Cross-validation**: 3-fold CV for robust evaluation
- **Feature scaling**: StandardScaler for regulatory factors
- **Data integration**: Averaging across molecular features per sample
- **Sample matching**: Automatic sample name standardization
- **Memory optimization**: Efficient handling of large datasets
- **Time series support**: Analysis across multiple time points and conditions
- **Network analysis**: Regulatory relationship density and pattern analysis

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
6. **Network Analysis**: Construct and analyze regulatory networks
7. **Time Series Analysis**: Examine temporal dynamics across conditions
8. **Visualization**: Generate comprehensive plots and summaries
9. **Results Export**: Save analysis results and visualizations

## 🎯 Key Findings

Based on the comprehensive analysis:

### Regulatory Hierarchy
1. **miRNAs**: Primary regulators with strong, direct effects (r = 0.170)
2. **lncRNAs**: Complex modulators with widespread influence (40,652 relationships)
3. **DNA Methylation**: Epigenetic controllers with subtle, long-term effects

### Network Characteristics
- **miRNA-Gene**: High regulatory density with repressive patterns
- **lncRNA-Gene**: Most extensive network with diverse mechanisms
- **Methylation-Gene**: Moderate density with context-dependent regulation

### Temporal Dynamics
- **4 time points analyzed** across **10 conditions**
- **Dynamic regulation patterns** identified
- **Condition-specific temporal responses** observed

## 🤝 Contributing

This pipeline is designed for bioinformatics research and can be extended with:
- Additional regulatory layers (e.g., transcription factors, chromatin accessibility)
- New machine learning models (e.g., deep learning, graph neural networks)
- Enhanced visualization options and interactive plots
- Integration with other omics data types (e.g., proteomics, metabolomics)
- Time series analysis improvements and temporal modeling
- Advanced network analysis and community detection algorithms

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
2. Review the `output/regulation_analysis/COMPREHENSIVE_ANALYSIS_SUMMARY_REPORT.md` for detailed results
3. Examine the output files for analysis results
4. Review the code comments for implementation details