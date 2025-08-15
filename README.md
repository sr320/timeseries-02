# Timeseries-02: Multi-Omics Gene Regulation Analysis Pipeline

A comprehensive bioinformatics pipeline for analyzing the influence of miRNA, lncRNA, and DNA methylation on gene expression using multi-omics data from 40 samples.

## 🧬 Project Overview

This repository contains a sophisticated analysis pipeline that examines how three regulatory layers influence gene expression patterns:
- **miRNA Expression** - microRNA regulatory molecules
- **lncRNA Expression** - long non-coding RNA molecules  
- **DNA Methylation** - epigenetic modifications at CpG sites

The pipeline applies multiple machine learning approaches to quantify regulatory influence and identify key regulatory relationships.

## 📁 Repository Structure

```
timeseries-02/
├── code/                          # Analysis scripts and pipeline
│   ├── analyze_gene_regulation.py    # Comprehensive regression analysis
│   ├── simple_correlation_analysis.py # Correlation-based analysis
│   ├── test_data_loading.py          # Data loading tests
│   ├── run_analysis.sh               # Automated pipeline runner
│   ├── requirements.txt              # Python dependencies
│   └── README.md                     # Detailed code documentation
├── data/                          # Input data files
│   ├── apul-gene_count_matrix.csv    # Gene expression (RNA-seq)
│   ├── Apul_miRNA_counts_formatted.txt    # miRNA expression
│   ├── Apul_lncRNA_counts_filtered.txt    # lncRNA expression
│   └── merged-WGBS-CpG-counts_filtered.csv # DNA methylation
├── output/                        # Analysis results and visualizations
└── README.md                      # This file
```

## 🔬 Analysis Methods

### 1. Data Preprocessing
- Multi-format data loading and cleaning
- Sample name standardization across datasets
- Log2 transformation for count data
- Quality filtering and zero-expression removal

### 2. Gene Selection
- Top 20 most variable genes selected for analysis
- Focus on genes likely to show strong regulatory patterns

### 3. Regulatory Influence Analysis
- **Multiple Linear Regression**: Basic linear relationships
- **Ridge Regression**: Regularized linear regression
- **Lasso Regression**: Sparse regression with feature selection
- **Random Forest**: Non-linear ensemble method
- **Correlation Analysis**: Simple correlation-based approach

### 4. Cross-Validation
- 3-fold cross-validation for robust model evaluation
- Feature scaling and standardization
- Comprehensive statistical evaluation

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
python simple_correlation_analysis.py
python analyze_gene_regulation.py
```

## 📊 Output Files

The pipeline generates comprehensive results including:
- **Statistical Results**: R² scores for each gene and model
- **Visualizations**: Heatmaps, distributions, and comparisons
- **Integrated Dataset**: Combined multi-omics dataset
- **Analysis Summary**: Top-performing genes and models

## 🧪 Data Sources

- **Gene Expression**: RNA-seq count matrix from 40 samples
- **miRNA**: Small RNA sequencing data
- **lncRNA**: Long non-coding RNA expression profiles
- **DNA Methylation**: Whole genome bisulfite sequencing (WGBS)

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

## 🛠️ Technical Details

- **Cross-validation**: 3-fold CV for robust evaluation
- **Feature scaling**: StandardScaler for regulatory factors
- **Data integration**: Averaging across molecular features per sample
- **Sample matching**: Automatic sample name standardization
- **Memory optimization**: Efficient handling of large datasets

## 📚 Dependencies

- **pandas** ≥1.3.0: Data manipulation and analysis
- **numpy** ≥1.21.0: Numerical computing
- **matplotlib** ≥3.4.0: Plotting and visualization
- **seaborn** ≥0.11.0: Statistical data visualization
- **scipy** ≥1.7.0: Scientific computing
- **scikit-learn** ≥1.0.0: Machine learning algorithms

## 🤝 Contributing

This pipeline is designed for bioinformatics research and can be extended with:
- Additional regulatory layers (e.g., transcription factors)
- New machine learning models
- Enhanced visualization options
- Integration with other omics data types

## 📄 License

[Add your license information here]

## 👥 Authors

[Add author information here]

## 📖 Citation

If you use this pipeline in your research, please cite:
[Add citation information here]