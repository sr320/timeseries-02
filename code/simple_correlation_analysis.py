#!/usr/bin/env python3
"""
Simple correlation analysis of miRNA, lncRNA, and DNA methylation influence on gene expression.
This script provides a straightforward correlation-based approach to understand regulatory relationships.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_data():
    """Load all data matrices."""
    print("Loading data...")
    
    # Load cleaned data files
    gene_data = pd.read_csv('../data/gene_counts_cleaned.csv')
    mirna_data = pd.read_csv('../data/mirna_counts_cleaned.csv')
    lncrna_data = pd.read_csv('../data/lncrna_counts_cleaned.csv')
    methylation_data = pd.read_csv('../data/methylation_counts_cleaned.csv')
    
    return gene_data, mirna_data, lncrna_data, methylation_data

def clean_gene_data(gene_data):
    """Clean gene expression data."""
    print("Processing gene expression data...")
    
    # Set gene_id as index
    gene_data.set_index('gene_id', inplace=True)
    
    # Log2 transform (add 1 to avoid log(0))
    gene_data_log = np.log2(gene_data + 1)
    
    return gene_data, gene_data_log

def clean_mirna_data(mirna_data):
    """Clean miRNA data."""
    print("Processing miRNA data...")
    
    # Set miRNA name as index
    mirna_data.set_index('Name', inplace=True)
    
    # Log2 transform
    mirna_data_log = np.log2(mirna_data + 1)
    
    return mirna_data, mirna_data_log

def clean_lncrna_data(lncrna_data):
    """Clean lncRNA data."""
    print("Processing lncRNA data...")
    
    # Set Geneid as index
    lncrna_data.set_index('Geneid', inplace=True)
    
    # Log2 transform
    lncrna_counts_log = np.log2(lncrna_data + 1)
    
    return lncrna_data, lncrna_counts_log

def clean_methylation_data(methylation_data):
    """Clean DNA methylation data."""
    print("Processing DNA methylation data...")
    
    # Set CpG as index
    methylation_data.set_index('CpG', inplace=True)
    
    return methylation_data

def select_genes(gene_data_log, n_genes=20):
    """Select top genes by variance."""
    print(f"Selecting top {n_genes} genes by variance...")
    
    gene_variance = gene_data_log.var(axis=1)
    top_genes = gene_variance.nlargest(n_genes).index.tolist()
    
    print(f"Selected genes: {top_genes[:5]}... (showing first 5)")
    
    return top_genes

def create_correlation_matrix(gene_data_log, mirna_data_log, lncrna_data_log, 
                            methylation_data, selected_genes, sample_mapping):
    """Create correlation matrix between genes and regulatory factors."""
    print("Creating correlation matrix...")
    
    # Get available samples
    available_samples = list(sample_mapping.keys())
    
    # Initialize results
    correlation_results = {}
    
    for gene in selected_genes:
        if gene in gene_data_log.index:
            print(f"Analyzing correlations for {gene}...")
            
            # Get gene expression for available samples
            gene_expr = {}
            for sample in available_samples:
                if sample in gene_data_log.columns:
                    gene_expr[sample] = gene_data_log.loc[gene, sample]
            
            if len(gene_expr) > 0:
                # Calculate correlations with regulatory factors
                correlations = {}
                
                # miRNA correlation (average across all miRNAs)
                mirna_corrs = []
                for sample in gene_expr.keys():
                    if sample in mirna_data_log.columns:
                        mirna_avg = mirna_data_log[sample].mean()
                        if not np.isnan(mirna_avg) and not np.isnan(gene_expr[sample]):
                            mirna_corrs.append((mirna_avg, gene_expr[sample]))
                
                if len(mirna_corrs) > 2:
                    mirna_values, gene_values = zip(*mirna_corrs)
                    mirna_corr, mirna_p = stats.pearsonr(mirna_values, gene_values)
                    correlations['miRNA'] = {'correlation': mirna_corr, 'p_value': mirna_p}
                
                # lncRNA correlation (average across all lncRNAs)
                lncrna_corrs = []
                for sample in gene_expr.keys():
                    if sample in lncrna_data_log.columns:
                        lncrna_avg = lncrna_data_log[sample].mean()
                        if not np.isnan(lncrna_avg) and not np.isnan(gene_expr[sample]):
                            lncrna_corrs.append((lncrna_avg, gene_expr[sample]))
                
                if len(lncrna_corrs) > 2:
                    lncrna_values, gene_values = zip(*lncrna_corrs)
                    lncrna_corr, lncrna_p = stats.pearsonr(lncrna_values, gene_values)
                    correlations['lncRNA'] = {'correlation': lncrna_corr, 'p_value': lncrna_p}
                
                # DNA methylation correlation (average across all CpGs)
                methylation_corrs = []
                for sample in gene_expr.keys():
                    if sample in methylation_data.columns:
                        methylation_avg = methylation_data[sample].mean()
                        if not np.isnan(methylation_avg) and not np.isnan(gene_expr[sample]):
                            methylation_corrs.append((methylation_avg, gene_expr[sample]))
                
                if len(methylation_corrs) > 2:
                    methylation_values, gene_values = zip(*methylation_corrs)
                    methylation_corr, methylation_p = stats.pearsonr(methylation_values, gene_values)
                    correlations['DNA_methylation'] = {'correlation': methylation_corr, 'p_value': methylation_p}
                
                correlation_results[gene] = correlations
    
    return correlation_results

def create_visualizations(correlation_results, output_dir):
    """Create correlation visualizations."""
    print("Creating visualizations...")
    
    # 1. Correlation heatmap
    plt.figure(figsize=(12, 8))
    
    # Prepare data for heatmap
    genes = list(correlation_results.keys())
    regulatory_factors = ['miRNA', 'lncRNA', 'DNA_methylation']
    
    heatmap_data = []
    for gene in genes:
        row = []
        for factor in regulatory_factors:
            if factor in correlation_results[gene]:
                row.append(correlation_results[gene][factor]['correlation'])
            else:
                row.append(np.nan)
        heatmap_data.append(row)
    
    heatmap_df = pd.DataFrame(heatmap_data, 
                            index=genes,
                            columns=regulatory_factors)
    
    sns.heatmap(heatmap_df, annot=True, cmap='RdBu_r', center=0,
               fmt='.3f', cbar_kws={'label': 'Correlation Coefficient'})
    plt.title('Gene-Regulatory Factor Correlations')
    plt.xlabel('Regulatory Factor')
    plt.ylabel('Gene')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/correlation_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Correlation distribution
    plt.figure(figsize=(15, 5))
    
    for i, factor in enumerate(regulatory_factors):
        plt.subplot(1, 3, i+1)
        
        correlations = []
        for gene in genes:
            if factor in correlation_results[gene]:
                correlations.append(correlation_results[gene][factor]['correlation'])
        
        if correlations:
            plt.hist(correlations, bins=15, alpha=0.7, edgecolor='black')
            plt.xlabel('Correlation Coefficient')
            plt.ylabel('Number of Genes')
            plt.title(f'{factor} Correlations')
            plt.axvline(np.mean(correlations), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(correlations):.3f}')
            plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/correlation_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Significance analysis
    plt.figure(figsize=(10, 6))
    
    significant_corrs = []
    factor_names = []
    
    for factor in regulatory_factors:
        for gene in genes:
            if factor in correlation_results[gene]:
                corr = correlation_results[gene][factor]['correlation']
                p_val = correlation_results[gene][factor]['p_value']
                if p_val < 0.05:  # Significant correlation
                    significant_corrs.append(corr)
                    factor_names.append(factor)
    
    if significant_corrs:
        plt.hist(significant_corrs, bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Correlation Coefficient (Significant only, p < 0.05)')
        plt.ylabel('Number of Correlations')
        plt.title('Distribution of Significant Correlations')
        plt.axvline(np.mean(significant_corrs), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(significant_corrs):.3f}')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{output_dir}/significant_correlations.png', dpi=300, bbox_inches='tight')
        plt.close()

def save_results(correlation_results, output_dir):
    """Save correlation results."""
    print("Saving results...")
    
    # Create results DataFrame
    results_data = []
    
    for gene in correlation_results:
        for factor, data in correlation_results[gene].items():
            results_data.append({
                'Gene': gene,
                'Regulatory_Factor': factor,
                'Correlation': data['correlation'],
                'P_Value': data['p_value'],
                'Significant': data['p_value'] < 0.05
            })
    
    results_df = pd.DataFrame(results_data)
    results_df.to_csv(f'{output_dir}/correlation_results.csv', index=False)
    
    # Create summary
    with open(f'{output_dir}/correlation_summary.txt', 'w') as f:
        f.write("CORRELATION ANALYSIS SUMMARY\n")
        f.write("=" * 40 + "\n\n")
        
        f.write(f"Total genes analyzed: {len(correlation_results)}\n")
        f.write(f"Total correlations calculated: {len(results_data)}\n")
        f.write(f"Significant correlations (p < 0.05): {len(results_df[results_df['Significant'] == True])}\n\n")
        
        # Summary by regulatory factor
        for factor in ['miRNA', 'lncRNA', 'DNA_methylation']:
            factor_data = results_df[results_df['Regulatory_Factor'] == factor]
            if len(factor_data) > 0:
                f.write(f"{factor}:\n")
                f.write(f"  Mean correlation: {factor_data['Correlation'].mean():.3f}\n")
                f.write(f"  Significant correlations: {len(factor_data[factor_data['Significant'] == True])}\n")
                f.write(f"  Strong correlations (|r| > 0.5): {len(factor_data[abs(factor_data['Correlation']) > 0.5])}\n\n")
        
        # Top correlations
        f.write("TOP 5 STRONGEST CORRELATIONS:\n")
        f.write("-" * 30 + "\n")
        top_corrs = results_df.nlargest(5, 'Correlation')
        for i, row in top_corrs.iterrows():
            f.write(f"{row['Gene']} - {row['Regulatory_Factor']}: r = {row['Correlation']:.3f} (p = {row['P_Value']:.3e})\n")
    
    print(f"Results saved to {output_dir}/")

def main():
    """Main analysis pipeline."""
    print("Starting Simple Correlation Analysis")
    print("=" * 40)
    
    # Create output directory
    import os
    output_dir = '../output'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Load data
        gene_data, mirna_data, lncrna_data, methylation_data = load_data()
        
        # Process data
        gene_data, gene_data_log = clean_gene_data(gene_data)
        mirna_data, mirna_data_log = clean_mirna_data(mirna_data)
        lncrna_data, lncrna_data_log = clean_lncrna_data(lncrna_data)
        methylation_data = clean_methylation_data(methylation_data)
        
        # Get available samples (all should be the same now)
        available_samples = list(gene_data.columns)
        print(f"Available samples: {len(available_samples)}")
        
        # Select genes
        selected_genes = select_genes(gene_data_log, n_genes=20)
        
        # Calculate correlations
        correlation_results = create_correlation_matrix(
            gene_data_log, mirna_data_log, lncrna_data_log,
            methylation_data, selected_genes, {sample: sample for sample in available_samples}
        )
        
        # Create visualizations
        create_visualizations(correlation_results, output_dir)
        
        # Save results
        save_results(correlation_results, output_dir)
        
        print("\nCorrelation analysis completed successfully!")
        print(f"Results saved to {output_dir}/")
        
    except Exception as e:
        print(f"Error in analysis pipeline: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
