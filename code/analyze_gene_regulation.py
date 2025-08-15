#!/usr/bin/env python3
"""
Analysis of miRNA, lncRNA, and DNA methylation influence on gene expression.
This script analyzes the regulatory effects of multiple molecular layers on gene expression.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def load_and_preprocess_data():
    """Load and preprocess all data matrices."""
    print("Loading data matrices...")
    
    # Load gene expression data
    print("Loading gene expression data...")
    gene_data = pd.read_csv('data/apul-gene_count_matrix.csv')
    
    # Load miRNA data
    print("Loading miRNA data...")
    mirna_data = pd.read_csv('data/Apul_miRNA_counts_formatted.txt', sep='\t')
    
    # Load lncRNA data
    print("Loading lncRNA data...")
    lncrna_data = pd.read_csv('data/Apul_lncRNA_counts_filtered.txt', sep='\t', comment='#')
    
    # Load DNA methylation data
    print("Loading DNA methylation data...")
    methylation_data = pd.read_csv('data/merged-WGBS-CpG-counts_filtered.csv')
    
    return gene_data, mirna_data, lncrna_data, methylation_data

def clean_gene_data(gene_data):
    """Clean and prepare gene expression data."""
    print("Cleaning gene expression data...")
    
    # Remove duplicate columns (some samples appear twice)
    gene_data = gene_data.loc[:, ~gene_data.columns.duplicated()]
    
    # Set gene_id as index
    gene_data.set_index('gene_id', inplace=True)
    
    # Remove genes with all zero counts
    gene_data = gene_data.loc[(gene_data != 0).any(axis=1)]
    
    # Log2 transform (add 1 to avoid log(0))
    gene_data_log = np.log2(gene_data + 1)
    
    return gene_data, gene_data_log

def clean_mirna_data(mirna_data):
    """Clean and prepare miRNA data."""
    print("Cleaning miRNA data...")
    
    # Set miRNA name as index
    mirna_data.set_index('Name', inplace=True)
    
    # Remove rows with all zero counts
    mirna_data = mirna_data.loc[(mirna_data != 0).any(axis=1)]
    
    # Log2 transform
    mirna_data_log = np.log2(mirna_data + 1)
    
    return mirna_data, mirna_data_log

def clean_lncrna_data(lncrna_data):
    """Clean and prepare lncRNA data."""
    print("Cleaning lncRNA data...")
    
    # Remove metadata columns and keep only count columns
    count_cols = [col for col in lncrna_data.columns if 'ACR-' in col]
    lncrna_counts = lncrna_data[['Geneid'] + count_cols].copy()
    
    # Set Geneid as index
    lncrna_counts.set_index('Geneid', inplace=True)
    
    # Remove rows with all zero counts
    lncrna_counts = lncrna_counts.loc[(lncrna_counts != 0).any(axis=1)]
    
    # Log2 transform
    lncrna_counts_log = np.log2(lncrna_counts + 1)
    
    return lncrna_counts, lncrna_counts_log

def clean_methylation_data(methylation_data):
    """Clean and prepare DNA methylation data."""
    print("Cleaning DNA methylation data...")
    
    # Set CpG as index
    methylation_data.set_index('CpG', inplace=True)
    
    # Remove rows with all zero values
    methylation_data = methylation_data.loc[(methylation_data != 0).any(axis=1)]
    
    # Fill NaN values with 0
    methylation_data = methylation_data.fillna(0)
    
    return methylation_data

def standardize_sample_names(data_dict):
    """Standardize sample names across all datasets."""
    print("Standardizing sample names...")
    
    # Extract sample names from gene data (most complete)
    gene_samples = [col for col in data_dict['gene'].columns if 'ACR-' in col]
    
    # Standardize sample names
    sample_mapping = {}
    for sample in gene_samples:
        # Extract ACR-XXX-TPX pattern
        if 'ACR-' in sample and 'TP' in sample:
            parts = sample.split('-')
            if len(parts) >= 3:
                acr_num = parts[1]
                tp_part = parts[2]
                if tp_part.startswith('TP'):
                    tp_num = tp_part[2:]
                    standardized = f"ACR-{acr_num}-TP{tp_num}"
                    sample_mapping[sample] = standardized
    
    return sample_mapping

def select_subset_genes(gene_data_log, n_genes=20):
    """Select a subset of genes for analysis."""
    print(f"Selecting {n_genes} genes for analysis...")
    
    # Select genes with highest variance (most variable)
    gene_variance = gene_data_log.var(axis=1)
    top_genes = gene_variance.nlargest(n_genes).index.tolist()
    
    print(f"Selected genes: {top_genes[:5]}... (showing first 5)")
    
    return top_genes

def create_integrated_dataset(gene_data_log, mirna_data_log, lncrna_data_log, 
                            methylation_data, sample_mapping, selected_genes):
    """Create integrated dataset with all molecular layers."""
    print("Creating integrated dataset...")
    
    # Get standardized sample names
    std_samples = list(set(sample_mapping.values()))
    std_samples.sort()
    
    # Initialize integrated dataset
    integrated_data = {}
    
    # Add gene expression data
    for gene in selected_genes:
        if gene in gene_data_log.index:
            gene_expr = gene_data_log.loc[gene]
            # Map to standardized sample names
            mapped_expr = {}
            for orig_sample, std_sample in sample_mapping.items():
                if orig_sample in gene_expr.index:
                    mapped_expr[std_sample] = gene_expr[orig_sample]
            
            integrated_data[f"gene_{gene}"] = mapped_expr
    
    # Add miRNA data (average across all miRNAs per sample)
    mirna_avg = mirna_data_log.mean(axis=0)
    for sample in std_samples:
        if sample in mirna_avg.index:
            integrated_data[f"mirna_avg"] = mirna_avg.to_dict()
    
    # Add lncRNA data (average across all lncRNAs per sample)
    lncrna_avg = lncrna_data_log.mean(axis=0)
    for sample in std_samples:
        if sample in lncrna_avg.index:
            integrated_data[f"lncrna_avg"] = lncrna_avg.to_dict()
    
    # Add methylation data (average across all CpGs per sample)
    methylation_avg = methylation_data.mean(axis=0)
    for sample in std_samples:
        if sample in methylation_avg.index:
            integrated_data[f"methylation_avg"] = methylation_avg.to_dict()
    
    # Convert to DataFrame
    integrated_df = pd.DataFrame(integrated_data)
    
    # Transpose to have samples as rows
    integrated_df = integrated_df.T
    
    return integrated_df

def analyze_regulatory_influence(integrated_df):
    """Analyze the influence of regulatory factors on gene expression."""
    print("Analyzing regulatory influence...")
    
    # Separate target variables (genes) and features (regulatory factors)
    gene_cols = [col for col in integrated_df.columns if col.startswith('gene_')]
    feature_cols = [col for col in integrated_df.columns if not col.startswith('gene_')]
    
    results = {}
    
    for gene_col in gene_cols:
        print(f"Analyzing {gene_col}...")
        
        # Prepare data for this gene
        y = integrated_df.loc[gene_col].dropna()
        X = integrated_df.loc[feature_cols, y.index].T
        
        if len(y) > 0 and len(X) > 0:
            # Standardize features
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            
            # Multiple linear regression
            try:
                lr_model = LinearRegression()
                lr_scores = cross_val_score(lr_model, X_scaled, y, cv=3, scoring='r2')
                lr_r2 = np.mean(lr_scores)
                
                # Ridge regression
                ridge_model = Ridge(alpha=1.0)
                ridge_scores = cross_val_score(ridge_model, X_scaled, y, cv=3, scoring='r2')
                ridge_r2 = np.mean(ridge_scores)
                
                # Lasso regression
                lasso_model = Lasso(alpha=0.1)
                lasso_scores = cross_val_score(lasso_model, X_scaled, y, cv=3, scoring='r2')
                lasso_r2 = np.mean(lasso_scores)
                
                # Random Forest
                rf_model = RandomForestRegressor(n_estimators=100, random_state=42)
                rf_scores = cross_val_score(rf_model, X_scaled, y, cv=3, scoring='r2')
                rf_r2 = np.mean(rf_scores)
                
                results[gene_col] = {
                    'linear_r2': lr_r2,
                    'ridge_r2': ridge_r2,
                    'lasso_r2': lasso_r2,
                    'random_forest_r2': rf_r2,
                    'mean_r2': np.mean([lr_r2, ridge_r2, lasso_r2, rf_r2])
                }
                
            except Exception as e:
                print(f"Error analyzing {gene_col}: {e}")
                results[gene_col] = {
                    'linear_r2': np.nan,
                    'ridge_r2': np.nan,
                    'lasso_r2': np.nan,
                    'random_forest_r2': np.nan,
                    'mean_r2': np.nan
                }
    
    return results

def create_visualizations(integrated_df, results, output_dir):
    """Create comprehensive visualizations."""
    print("Creating visualizations...")
    
    # 1. Overall regulatory influence heatmap
    plt.figure(figsize=(12, 8))
    
    # Prepare data for heatmap
    gene_cols = [col for col in integrated_df.index if col.startswith('gene_')]
    feature_cols = [col for col in integrated_df.index if not col.startswith('gene_')]
    
    heatmap_data = []
    for gene in gene_cols:
        if gene in results:
            row = [
                results[gene]['linear_r2'],
                results[gene]['ridge_r2'],
                results[gene]['lasso_r2'],
                results[gene]['random_forest_r2']
            ]
            heatmap_data.append(row)
    
    if heatmap_data:
        heatmap_df = pd.DataFrame(heatmap_data, 
                                index=gene_cols,
                                columns=['Linear', 'Ridge', 'Lasso', 'Random Forest'])
        
        sns.heatmap(heatmap_df, annot=True, cmap='RdYlBu_r', center=0,
                   fmt='.3f', cbar_kws={'label': 'R² Score'})
        plt.title('Regulatory Influence Analysis: Model Performance by Gene')
        plt.xlabel('Regression Model')
        plt.ylabel('Gene')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/regulatory_influence_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 2. Distribution of R² scores
    plt.figure(figsize=(10, 6))
    
    r2_scores = []
    for gene in results:
        r2_scores.append(results[gene]['mean_r2'])
    
    plt.hist(r2_scores, bins=20, alpha=0.7, edgecolor='black')
    plt.xlabel('Mean R² Score')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of Regulatory Influence Scores')
    plt.axvline(np.mean(r2_scores), color='red', linestyle='--', 
                label=f'Mean: {np.mean(r2_scores):.3f}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/r2_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Model comparison
    plt.figure(figsize=(10, 6))
    
    model_scores = {
        'Linear': [results[gene]['linear_r2'] for gene in results if not np.isnan(results[gene]['linear_r2'])],
        'Ridge': [results[gene]['ridge_r2'] for gene in results if not np.isnan(results[gene]['ridge_r2'])],
        'Lasso': [results[gene]['lasso_r2'] for gene in results if not np.isnan(results[gene]['lasso_r2'])],
        'Random Forest': [results[gene]['random_forest_r2'] for gene in results if not np.isnan(results[gene]['random_forest_r2'])]
    }
    
    plt.boxplot(model_scores.values(), labels=model_scores.keys())
    plt.ylabel('R² Score')
    plt.title('Model Performance Comparison')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/model_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Sample correlation heatmap
    plt.figure(figsize=(12, 10))
    
    # Get sample correlation matrix
    sample_corr = integrated_df.T.corr()
    
    sns.heatmap(sample_corr, cmap='coolwarm', center=0, 
               square=True, cbar_kws={'label': 'Correlation'})
    plt.title('Sample Correlation Matrix')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/sample_correlation.png', dpi=300, bbox_inches='tight')
    plt.close()

def save_results(results, integrated_df, output_dir):
    """Save analysis results to files."""
    print("Saving results...")
    
    # Save results summary
    results_df = pd.DataFrame(results).T
    results_df.to_csv(f'{output_dir}/regulatory_influence_results.csv')
    
    # Save integrated dataset
    integrated_df.to_csv(f'{output_dir}/integrated_dataset.csv')
    
    # Create summary statistics
    summary_stats = {
        'total_genes_analyzed': len(results),
        'mean_linear_r2': np.mean([results[gene]['linear_r2'] for gene in results if not np.isnan(results[gene]['linear_r2'])]),
        'mean_ridge_r2': np.mean([results[gene]['ridge_r2'] for gene in results if not np.isnan(results[gene]['ridge_r2'])]),
        'mean_lasso_r2': np.mean([results[gene]['lasso_r2'] for gene in results if not np.isnan(results[gene]['lasso_r2'])]),
        'mean_rf_r2': np.mean([results[gene]['random_forest_r2'] for gene in results if not np.isnan(results[gene]['random_forest_r2'])]),
        'best_performing_model': max(['Linear', 'Ridge', 'Lasso', 'Random Forest'], 
                                   key=lambda x: np.mean([results[gene][f'{x.lower().replace(" ", "_")}_r2'] 
                                                        for gene in results if not np.isnan(results[gene][f'{x.lower().replace(" ", "_")}_r2'])])),
        'genes_with_high_influence': len([gene for gene in results if results[gene]['mean_r2'] > 0.5])
    }
    
    # Save summary
    with open(f'{output_dir}/analysis_summary.txt', 'w') as f:
        f.write("GENE REGULATION ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total genes analyzed: {summary_stats['total_genes_analyzed']}\n")
        f.write(f"Mean Linear R²: {summary_stats['mean_linear_r2']:.3f}\n")
        f.write(f"Mean Ridge R²: {summary_stats['mean_ridge_r2']:.3f}\n")
        f.write(f"Mean Lasso R²: {summary_stats['mean_lasso_r2']:.3f}\n")
        f.write(f"Mean Random Forest R²: {summary_stats['mean_rf_r2']:.3f}\n")
        f.write(f"Best performing model: {summary_stats['best_performing_model']}\n")
        f.write(f"Genes with high regulatory influence (R² > 0.5): {summary_stats['genes_with_high_influence']}\n")
        
        f.write("\nTOP 5 GENES BY REGULATORY INFLUENCE:\n")
        f.write("-" * 40 + "\n")
        sorted_genes = sorted(results.items(), key=lambda x: x[1]['mean_r2'], reverse=True)
        for i, (gene, scores) in enumerate(sorted_genes[:5]):
            f.write(f"{i+1}. {gene}: Mean R² = {scores['mean_r2']:.3f}\n")
    
    print(f"Results saved to {output_dir}/")

def main():
    """Main analysis pipeline."""
    print("Starting Gene Regulation Analysis Pipeline")
    print("=" * 50)
    
    # Create output directory
    import os
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Load data
        gene_data, mirna_data, lncrna_data, methylation_data = load_and_preprocess_data()
        
        # Clean data
        gene_data, gene_data_log = clean_gene_data(gene_data)
        mirna_data, mirna_data_log = clean_mirna_data(mirna_data)
        lncrna_data, lncrna_data_log = clean_lncrna_data(lncrna_data)
        methylation_data = clean_methylation_data(methylation_data)
        
        # Standardize sample names
        sample_mapping = standardize_sample_names({
            'gene': gene_data,
            'mirna': mirna_data,
            'lncrna': lncrna_data,
            'methylation': methylation_data
        })
        
        # Select subset of genes
        selected_genes = select_subset_genes(gene_data_log, n_genes=20)
        
        # Create integrated dataset
        integrated_df = create_integrated_dataset(
            gene_data_log, mirna_data_log, lncrna_data_log, 
            methylation_data, sample_mapping, selected_genes
        )
        
        # Analyze regulatory influence
        results = analyze_regulatory_influence(integrated_df)
        
        # Create visualizations
        create_visualizations(integrated_df, results, output_dir)
        
        # Save results
        save_results(results, integrated_df, output_dir)
        
        print("\nAnalysis completed successfully!")
        print(f"Results saved to {output_dir}/")
        
    except Exception as e:
        print(f"Error in analysis pipeline: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
