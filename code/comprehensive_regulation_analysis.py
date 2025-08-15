#!/usr/bin/env python3
"""
Comprehensive Analysis: How miRNA, lncRNA, and DNA methylation influence gene expression

This script performs integrated multi-omics analysis to understand regulatory mechanisms:
1. Correlation analysis between different data types
2. Time series analysis across time points
3. Regulatory network inference
4. Statistical modeling of expression regulation
5. Visualization of regulatory relationships
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
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

class ComprehensiveRegulationAnalysis:
    def __init__(self, data_dir="data/cleaned_datasets"):
        """Initialize the analysis with cleaned datasets."""
        self.data_dir = data_dir
        self.datasets = {}
        self.results = {}
        self.load_datasets()
        
    def load_datasets(self):
        """Load all cleaned datasets."""
        print("Loading cleaned datasets...")
        
        # Load gene expression data
        self.datasets['gene'] = pd.read_csv(f"{self.data_dir}/gene_counts_cleaned.csv", index_col=0)
        print(f"Loaded gene expression: {self.datasets['gene'].shape}")
        
        # Load lncRNA expression data
        self.datasets['lncrna'] = pd.read_csv(f"{self.data_dir}/lncrna_counts_cleaned.csv", index_col=0)
        print(f"Loaded lncRNA expression: {self.datasets['lncrna'].shape}")
        
        # Load miRNA expression data
        self.datasets['mirna'] = pd.read_csv(f"{self.data_dir}/mirna_counts_cleaned.csv", index_col=0)
        print(f"Loaded miRNA expression: {self.datasets['mirna'].shape}")
        
        # Load DNA methylation data
        self.datasets['methylation'] = pd.read_csv(f"{self.data_dir}/wgbs_counts_cleaned.csv", index_col=0)
        print(f"Loaded DNA methylation: {self.datasets['methylation'].shape}")
        
        # Verify sample alignment
        self.verify_sample_alignment()
        
    def verify_sample_alignment(self):
        """Verify that all datasets have the same sample structure."""
        sample_sets = [set(df.columns) for df in self.datasets.values()]
        if not all(sample_sets[0] == s for s in sample_sets):
            raise ValueError("Sample IDs are not aligned across datasets!")
        
        self.samples = sorted(list(sample_sets[0]))
        self.n_samples = len(self.samples)
        print(f"✓ All datasets aligned with {self.n_samples} samples")
        
        # Extract time points and conditions
        self.time_points = sorted(list(set([s.split('-')[-1] for s in self.samples])))
        self.conditions = sorted(list(set([s.split('-')[1] for s in self.samples])))
        print(f"Time points: {self.time_points}")
        print(f"Conditions: {self.conditions}")
        
    def correlation_analysis(self):
        """Perform comprehensive correlation analysis between data types."""
        print("\n" + "="*60)
        print("CORRELATION ANALYSIS")
        print("="*60)
        
        # 1. Gene vs lncRNA correlation
        print("Analyzing gene-lncRNA correlations...")
        gene_lncrna_corr = self.calculate_cross_correlations(
            self.datasets['gene'], 
            self.datasets['lncrna'], 
            'gene_lncrna'
        )
        
        # 2. Gene vs miRNA correlation
        print("Analyzing gene-miRNA correlations...")
        gene_mirna_corr = self.calculate_cross_correlations(
            self.datasets['gene'], 
            self.datasets['mirna'], 
            'gene_mirna'
        )
        
        # 3. Gene vs methylation correlation
        print("Analyzing gene-methylation correlations...")
        gene_meth_corr = self.calculate_cross_correlations(
            self.datasets['gene'], 
            self.datasets['methylation'], 
            'gene_methylation'
        )
        
        # 4. lncRNA vs miRNA correlation
        print("Analyzing lncRNA-miRNA correlations...")
        lncrna_mirna_corr = self.calculate_cross_correlations(
            self.datasets['lncrna'], 
            self.datasets['mirna'], 
            'lncrna_mirna'
        )
        
        # Store results
        self.results['correlations'] = {
            'gene_lncrna': gene_lncrna_corr,
            'gene_mirna': gene_mirna_corr,
            'gene_methylation': gene_meth_corr,
            'lncrna_mirna': lncrna_mirna_corr
        }
        
        return self.results['correlations']
    
    def calculate_cross_correlations(self, data1, data2, name):
        """Calculate correlations between two datasets."""
        print(f"  Calculating {name} correlations...")
        
        # Sample correlations (across features for each sample)
        sample_corrs = []
        for sample in self.samples:
            if sample in data1.columns and sample in data2.columns:
                # Get expression values for this sample from both datasets
                expr1 = data1[sample].values
                expr2 = data2[sample].values
                
                # Ensure both arrays have the same length by truncating to shorter one
                min_length = min(len(expr1), len(expr2))
                if min_length > 1:  # Need at least 2 points for correlation
                    corr, pval = pearsonr(expr1[:min_length], expr2[:min_length])
                    sample_corrs.append({'sample': sample, 'correlation': corr, 'p_value': pval})
        
        sample_corr_df = pd.DataFrame(sample_corrs)
        
        # Feature correlations (across samples for each feature)
        # For computational efficiency, sample a subset of features
        n_features1 = min(1000, len(data1))
        n_features2 = min(1000, len(data2))
        
        sampled_features1 = np.random.choice(data1.index, n_features1, replace=False)
        sampled_features2 = np.random.choice(data2.index, n_features2, replace=False)
        
        feature_corrs = []
        for feat1 in sampled_features1[:100]:  # Limit for computational efficiency
            for feat2 in sampled_features2[:100]:
                if feat1 in data1.index and feat2 in data2.index:
                    # Get expression values across samples
                    expr1 = data1.loc[feat1].values
                    expr2 = data2.loc[feat2].values
                    
                    # Ensure both arrays have the same length
                    min_length = min(len(expr1), len(expr2))
                    if min_length > 1:
                        corr, pval = pearsonr(expr1[:min_length], expr2[:min_length])
                        feature_corrs.append({
                            'feature1': feat1, 'feature2': feat2, 
                            'correlation': corr, 'p_value': pval
                        })
        
        feature_corr_df = pd.DataFrame(feature_corrs)
        
        return {
            'sample_correlations': sample_corr_df,
            'feature_correlations': feature_corr_df
        }
    
    def time_series_analysis(self):
        """Analyze changes across time points."""
        print("\n" + "="*60)
        print("TIME SERIES ANALYSIS")
        print("="*60)
        
        time_results = {}
        
        for condition in self.conditions:
            print(f"Analyzing time series for condition {condition}...")
            
            # Get samples for this condition
            condition_samples = [s for s in self.samples if f"ACR-{condition}-" in s]
            condition_samples = sorted(condition_samples)
            
            if len(condition_samples) == 4:  # All time points available
                time_results[condition] = {}
                
                # Analyze each data type across time
                for data_type, data in self.datasets.items():
                    condition_data = data[condition_samples]
                    
                    # Calculate temporal trends
                    trends = self.calculate_temporal_trends(condition_data)
                    time_results[condition][data_type] = trends
        
        self.results['time_series'] = time_results
        return time_results
    
    def calculate_temporal_trends(self, data):
        """Calculate temporal trends for a dataset."""
        trends = {}
        
        # Linear trend across time points
        time_values = np.array([1, 2, 3, 4])  # TP1, TP2, TP3, TP4
        
        for feature in data.index:
            expression_values = data.loc[feature].values
            
            if len(expression_values) == 4:
                # Linear regression
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    time_values, expression_values
                )
                
                trends[feature] = {
                    'slope': slope,
                    'intercept': intercept,
                    'r_squared': r_value**2,
                    'p_value': p_value,
                    'trend_direction': 'increasing' if slope > 0 else 'decreasing',
                    'trend_strength': abs(slope)
                }
        
        return trends
    
    def regulatory_network_analysis(self):
        """Infer regulatory networks between different data types."""
        print("\n" + "="*60)
        print("REGULATORY NETWORK ANALYSIS")
        print("="*60)
        
        # 1. miRNA-gene regulatory network
        print("Inferring miRNA-gene regulatory network...")
        mirna_gene_network = self.infer_mirna_gene_regulation()
        
        # 2. lncRNA-gene regulatory network
        print("Inferring lncRNA-gene regulatory network...")
        lncrna_gene_network = self.infer_lncrna_gene_regulation()
        
        # 3. Methylation-gene regulatory network
        print("Inferring methylation-gene regulatory network...")
        methylation_gene_network = self.infer_methylation_gene_regulation()
        
        # Store results
        self.results['regulatory_networks'] = {
            'mirna_gene': mirna_gene_network,
            'lncrna_gene': lncrna_gene_network,
            'methylation_gene': methylation_gene_network
        }
        
        return self.results['regulatory_networks']
    
    def infer_mirna_gene_regulation(self):
        """Infer miRNA-gene regulatory relationships."""
        # For miRNA-gene regulation, we typically expect negative correlations
        # as miRNAs often repress gene expression
        
        mirna_gene_network = []
        
        # Sample a subset for computational efficiency
        n_genes = min(500, len(self.datasets['gene']))
        n_mirnas = min(51, len(self.datasets['mirna']))
        
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        sampled_mirnas = np.random.choice(self.datasets['mirna'].index, n_mirnas, replace=False)
        
        for gene in sampled_genes:
            for mirna in sampled_mirnas:
                corr, pval = pearsonr(
                    self.datasets['gene'].loc[gene], 
                    self.datasets['mirna'].loc[mirna]
                )
                
                if pval < 0.05:  # Significant correlation
                    regulation_type = 'repression' if corr < 0 else 'activation'
                    mirna_gene_network.append({
                        'mirna': mirna,
                        'gene': gene,
                        'correlation': corr,
                        'p_value': pval,
                        'regulation_type': regulation_type,
                        'strength': abs(corr)
                    })
        
        return pd.DataFrame(mirna_gene_network)
    
    def infer_lncrna_gene_regulation(self):
        """Infer lncRNA-gene regulatory relationships."""
        lncrna_gene_network = []
        
        # Sample a subset for computational efficiency
        n_genes = min(500, len(self.datasets['gene']))
        n_lncrnas = min(500, len(self.datasets['lncrna']))
        
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        sampled_lncrnas = np.random.choice(self.datasets['lncrna'].index, n_lncrnas, replace=False)
        
        for gene in sampled_genes:
            for lncrna in sampled_lncrnas:
                corr, pval = pearsonr(
                    self.datasets['gene'].loc[gene], 
                    self.datasets['lncrna'].loc[lncrna]
                )
                
                if pval < 0.05:  # Significant correlation
                    regulation_type = 'activation' if corr > 0 else 'repression'
                    lncrna_gene_network.append({
                        'lncrna': lncrna,
                        'gene': gene,
                        'correlation': corr,
                        'p_value': pval,
                        'regulation_type': regulation_type,
                        'strength': abs(corr)
                    })
        
        return pd.DataFrame(lncrna_gene_network)
    
    def infer_methylation_gene_regulation(self):
        """Infer methylation-gene regulatory relationships."""
        methylation_gene_network = []
        
        # Sample a subset for computational efficiency
        n_genes = min(500, len(self.datasets['gene']))
        n_cpgs = min(249, len(self.datasets['methylation']))
        
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        sampled_cpgs = np.random.choice(self.datasets['methylation'].index, n_cpgs, replace=False)
        
        for gene in sampled_genes:
            for cpg in sampled_cpgs:
                corr, pval = pearsonr(
                    self.datasets['gene'].loc[gene], 
                    self.datasets['methylation'].loc[cpg]
                )
                
                if pval < 0.05:  # Significant correlation
                    regulation_type = 'repression' if corr < 0 else 'activation'
                    methylation_gene_network.append({
                        'cpg': cpg,
                        'gene': gene,
                        'correlation': corr,
                        'p_value': pval,
                        'regulation_type': regulation_type,
                        'strength': abs(corr)
                    })
        
        return pd.DataFrame(methylation_gene_network)
    
    def predictive_modeling(self):
        """Build predictive models for gene expression regulation."""
        print("\n" + "="*60)
        print("PREDICTIVE MODELING")
        print("="*60)
        
        # Prepare integrated dataset
        print("Preparing integrated dataset for modeling...")
        integrated_data = self.prepare_integrated_dataset()
        
        # Build models
        models = {}
        
        # 1. Linear Regression
        print("Building Linear Regression model...")
        models['linear'] = self.build_linear_model(integrated_data)
        
        # 2. Ridge Regression
        print("Building Ridge Regression model...")
        models['ridge'] = self.build_ridge_model(integrated_data)
        
        # 3. Random Forest
        print("Building Random Forest model...")
        models['random_forest'] = self.build_random_forest_model(integrated_data)
        
        self.results['predictive_models'] = models
        return models
    
    def prepare_integrated_dataset(self):
        """Prepare integrated dataset for predictive modeling."""
        # For modeling, we'll focus on a subset of genes and their regulators
        n_genes = min(200, len(self.datasets['gene']))
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        
        # Create integrated dataset
        integrated_data = {}
        
        for gene in sampled_genes:
            gene_expression = self.datasets['gene'].loc[gene]
            
            # Get corresponding regulatory features
            regulators = {}
            
            # miRNA regulators (top correlations)
            mirna_corrs = []
            for mirna in self.datasets['mirna'].index:
                corr, pval = pearsonr(gene_expression, self.datasets['mirna'].loc[mirna])
                if pval < 0.1:  # Relaxed threshold for feature selection
                    mirna_corrs.append((mirna, corr, pval))
            
            # Sort by correlation strength and take top 5
            mirna_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
            for mirna, corr, pval in mirna_corrs[:5]:
                regulators[f"mirna_{mirna}"] = self.datasets['mirna'].loc[mirna]
            
            # lncRNA regulators (top correlations)
            lncrna_corrs = []
            for lncrna in self.datasets['lncrna'].index:
                corr, pval = pearsonr(gene_expression, self.datasets['lncrna'].loc[lncrna])
                if pval < 0.1:
                    lncrna_corrs.append((lncrna, corr, pval))
            
            lncrna_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
            for lncrna, corr, pval in lncrna_corrs[:5]:
                regulators[f"lncrna_{lncrna}"] = self.datasets['lncrna'].loc[lncrna]
            
            # Methylation regulators (top correlations)
            meth_corrs = []
            for cpg in self.datasets['methylation'].index:
                corr, pval = pearsonr(gene_expression, self.datasets['methylation'].loc[cpg])
                if pval < 0.1:
                    meth_corrs.append((cpg, corr, pval))
            
            meth_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
            for cpg, corr, pval in meth_corrs[:5]:
                regulators[f"methylation_{cpg}"] = self.datasets['methylation'].loc[cpg]
            
            # Create feature matrix
            if regulators:
                feature_matrix = pd.DataFrame(regulators)
                feature_matrix['gene_expression'] = gene_expression
                integrated_data[gene] = feature_matrix
        
        return integrated_data
    
    def build_linear_model(self, integrated_data):
        """Build linear regression model."""
        results = {}
        
        for gene, data in integrated_data.items():
            if len(data) > 10:  # Need sufficient data
                X = data.drop('gene_expression', axis=1)
                y = data['gene_expression']
                
                # Split data
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=0.3, random_state=42
                )
                
                # Train model
                model = LinearRegression()
                model.fit(X_train, y_train)
                
                # Predictions
                y_pred = model.predict(X_test)
                
                # Evaluate
                r2 = r2_score(y_test, y_pred)
                mse = mean_squared_error(y_test, y_pred)
                
                results[gene] = {
                    'model': model,
                    'r2_score': r2,
                    'mse': mse,
                    'feature_importance': dict(zip(X.columns, model.coef_))
                }
        
        return results
    
    def build_ridge_model(self, integrated_data):
        """Build ridge regression model."""
        results = {}
        
        for gene, data in integrated_data.items():
            if len(data) > 10:
                X = data.drop('gene_expression', axis=1)
                y = data['gene_expression']
                
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=0.3, random_state=42
                )
                
                model = Ridge(alpha=1.0)
                model.fit(X_train, y_train)
                
                y_pred = model.predict(X_test)
                r2 = r2_score(y_test, y_pred)
                mse = mean_squared_error(y_test, y_pred)
                
                results[gene] = {
                    'model': model,
                    'r2_score': r2,
                    'mse': mse,
                    'feature_importance': dict(zip(X.columns, model.coef_))
                }
        
        return results
    
    def build_random_forest_model(self, integrated_data):
        """Build random forest model."""
        results = {}
        
        for gene, data in integrated_data.items():
            if len(data) > 10:
                X = data.drop('gene_expression', axis=1)
                y = data['gene_expression']
                
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=0.3, random_state=42
                )
                
                model = RandomForestRegressor(n_estimators=100, random_state=42)
                model.fit(X_train, y_train)
                
                y_pred = model.predict(X_test)
                r2 = r2_score(y_test, y_pred)
                mse = mean_squared_error(y_test, y_pred)
                
                results[gene] = {
                    'model': model,
                    'r2_score': r2,
                    'mse': mse,
                    'feature_importance': dict(zip(X.columns, model.feature_importances_))
                }
        
        return results
    
    def generate_visualizations(self):
        """Generate comprehensive visualizations."""
        print("\n" + "="*60)
        print("GENERATING VISUALIZATIONS")
        print("="*60)
        
        # Create output directory
        import os
        os.makedirs("output/regulation_analysis", exist_ok=True)
        
        # 1. Correlation heatmaps
        self.plot_correlation_heatmaps()
        
        # 2. Time series plots
        self.plot_time_series()
        
        # 3. Regulatory network plots
        self.plot_regulatory_networks()
        
        # 4. Model performance plots
        self.plot_model_performance()
        
        print("✓ All visualizations saved to output/regulation_analysis/")
    
    def plot_correlation_heatmaps(self):
        """Plot correlation heatmaps."""
        # Sample correlations across samples
        sample_corrs = {}
        for data_type in ['gene', 'lncrna', 'mirna', 'methylation']:
            sample_corrs[data_type] = self.datasets[data_type].T.corr()
        
        # Plot sample correlation heatmap
        plt.figure(figsize=(12, 10))
        plt.subplot(2, 2, 1)
        sns.heatmap(sample_corrs['gene'].iloc[:20, :20], cmap='RdBu_r', center=0)
        plt.title('Gene Expression Sample Correlations')
        
        plt.subplot(2, 2, 2)
        sns.heatmap(sample_corrs['lncrna'].iloc[:20, :20], cmap='RdBu_r', center=0)
        plt.title('lncRNA Expression Sample Correlations')
        
        plt.subplot(2, 2, 3)
        sns.heatmap(sample_corrs['mirna'].iloc[:20, :20], cmap='RdBu_r', center=0)
        plt.title('miRNA Expression Sample Correlations')
        
        plt.subplot(2, 2, 4)
        sns.heatmap(sample_corrs['methylation'].iloc[:20, :20], cmap='RdBu_r', center=0)
        plt.title('Methylation Sample Correlations')
        
        plt.tight_layout()
        plt.savefig("output/regulation_analysis/sample_correlation_heatmaps.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_time_series(self):
        """Plot time series analysis."""
        # Plot average expression across time points for each condition
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.ravel()
        
        for i, data_type in enumerate(['gene', 'lncrna', 'mirna', 'methylation']):
            data = self.datasets[data_type]
            
            # Calculate mean expression per time point per condition
            time_means = {}
            for condition in self.conditions:
                condition_samples = [s for s in self.samples if f"ACR-{condition}-" in s]
                condition_samples = sorted(condition_samples)
                
                if len(condition_samples) == 4:
                    time_means[condition] = data[condition_samples].mean(axis=0)
            
            # Plot
            for condition, means in time_means.items():
                time_points = [1, 2, 3, 4]
                axes[i].plot(time_points, means.values, marker='o', label=f'ACR-{condition}')
            
            axes[i].set_title(f'{data_type.upper()} Expression Across Time')
            axes[i].set_xlabel('Time Point')
            axes[i].set_ylabel('Mean Expression')
            axes[i].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            axes[i].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig("output/regulation_analysis/time_series_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_regulatory_networks(self):
        """Plot regulatory network summaries."""
        if 'regulatory_networks' not in self.results:
            print("No regulatory networks to plot. Run regulatory_network_analysis() first.")
            return
        
        # Plot regulation type distributions
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        # miRNA-gene regulation
        if not self.results['regulatory_networks']['mirna_gene'].empty:
            regulation_counts = self.results['regulatory_networks']['mirna_gene']['regulation_type'].value_counts()
            axes[0].pie(regulation_counts.values, labels=regulation_counts.index, autopct='%1.1f%%')
            axes[0].set_title('miRNA-Gene Regulation Types')
        
        # lncRNA-gene regulation
        if not self.results['regulatory_networks']['lncrna_gene'].empty:
            regulation_counts = self.results['regulatory_networks']['lncrna_gene']['regulation_type'].value_counts()
            axes[1].pie(regulation_counts.values, labels=regulation_counts.index, autopct='%1.1f%%')
            axes[1].set_title('lncRNA-Gene Regulation Types')
        
        # Methylation-gene regulation
        if not self.results['regulatory_networks']['methylation_gene'].empty:
            regulation_counts = self.results['regulatory_networks']['methylation_gene']['regulation_type'].value_counts()
            axes[2].pie(regulation_counts.values, labels=regulation_counts.index, autopct='%1.1f%%')
            axes[2].set_title('Methylation-Gene Regulation Types')
        
        plt.tight_layout()
        plt.savefig("output/regulation_analysis/regulation_type_distributions.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_model_performance(self):
        """Plot model performance comparisons."""
        if 'predictive_models' not in self.results:
            print("No predictive models to plot. Run predictive_modeling() first.")
            return
        
        # Compare R² scores across models
        model_names = list(self.results['predictive_models'].keys())
        r2_scores = {model: [] for model in model_names}
        
        for model_name in model_names:
            for gene, results in self.results['predictive_models'][model_name].items():
                r2_scores[model_name].append(results['r2_score'])
        
        # Create box plot
        plt.figure(figsize=(10, 6))
        data_to_plot = [r2_scores[model] for model in model_names]
        plt.boxplot(data_to_plot, labels=model_names)
        plt.title('Model Performance Comparison (R² Scores)')
        plt.ylabel('R² Score')
        plt.grid(True, alpha=0.3)
        
        plt.savefig("output/regulation_analysis/model_performance_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    def save_results(self):
        """Save all analysis results."""
        print("\n" + "="*60)
        print("SAVING RESULTS")
        print("="*60)
        
        import os
        os.makedirs("output/regulation_analysis", exist_ok=True)
        
        # Save correlation results
        for corr_type, corr_data in self.results.get('correlations', {}).items():
            for data_type, data in corr_data.items():
                if isinstance(data, pd.DataFrame):
                    data.to_csv(f"output/regulation_analysis/{corr_type}_{data_type}.csv")
        
        # Save regulatory networks
        for network_type, network_data in self.results.get('regulatory_networks', {}).items():
            if isinstance(network_data, pd.DataFrame):
                network_data.to_csv(f"output/regulation_analysis/{network_type}_network.csv")
        
        # Save time series results
        if 'time_series' in self.results:
            import json
            # Convert numpy types to native Python types for JSON serialization
            time_series_serializable = {}
            for condition, data_types in self.results['time_series'].items():
                time_series_serializable[condition] = {}
                for data_type, trends in data_types.items():
                    time_series_serializable[condition][data_type] = {}
                    for feature, trend_data in trends.items():
                        time_series_serializable[condition][data_type][feature] = {
                            k: float(v) if isinstance(v, (np.float32, np.float64)) else v
                            for k, v in trend_data.items()
                        }
            
            with open("output/regulation_analysis/time_series_results.json", 'w') as f:
                json.dump(time_series_serializable, f, indent=2)
        
        # Save model results summary
        if 'predictive_models' in self.results:
            model_summary = {}
            for model_type, gene_results in self.results['predictive_models'].items():
                model_summary[model_type] = {
                    'mean_r2': np.mean([r['r2_score'] for r in gene_results.values()]),
                    'mean_mse': np.mean([r['mse'] for r in gene_results.values()]),
                    'n_genes': len(gene_results)
                }
            
            pd.DataFrame(model_summary).T.to_csv("output/regulation_analysis/model_performance_summary.csv")
        
        print("✓ All results saved to output/regulation_analysis/")
    
    def run_complete_analysis(self):
        """Run the complete comprehensive analysis."""
        print("="*80)
        print("COMPREHENSIVE REGULATION ANALYSIS")
        print("="*80)
        print("This analysis will determine how miRNA, lncRNA, and DNA methylation")
        print("influence gene expression using integrated multi-omics data.")
        print("="*80)
        
        # Run all analyses
        self.correlation_analysis()
        self.time_series_analysis()
        self.regulatory_network_analysis()
        self.predictive_modeling()
        
        # Generate visualizations
        self.generate_visualizations()
        
        # Save results
        self.save_results()
        
        # Print summary
        self.print_analysis_summary()
        
        print("\n" + "="*80)
        print("ANALYSIS COMPLETED SUCCESSFULLY!")
        print("="*80)
    
    def print_analysis_summary(self):
        """Print a comprehensive summary of the analysis results."""
        print("\n" + "="*60)
        print("ANALYSIS SUMMARY")
        print("="*60)
        
        # Correlation summary
        if 'correlations' in self.results:
            print("\nCORRELATION ANALYSIS RESULTS:")
            for corr_type, corr_data in self.results['correlations'].items():
                if 'sample_correlations' in corr_data:
                    sample_corrs = corr_data['sample_correlations']
                    mean_corr = sample_corrs['correlation'].mean()
                    print(f"  {corr_type}: Mean sample correlation = {mean_corr:.3f}")
        
        # Regulatory network summary
        if 'regulatory_networks' in self.results:
            print("\nREGULATORY NETWORK RESULTS:")
            for network_type, network_data in self.results['regulatory_networks'].items():
                if isinstance(network_data, pd.DataFrame) and not network_data.empty:
                    n_regulations = len(network_data)
                    print(f"  {network_type}: {n_regulations} regulatory relationships identified")
        
        # Model performance summary
        if 'predictive_models' in self.results:
            print("\nPREDICTIVE MODEL PERFORMANCE:")
            for model_type, gene_results in self.results['predictive_models'].items():
                if gene_results:
                    mean_r2 = np.mean([r['r2_score'] for r in gene_results.values()])
                    print(f"  {model_type}: Mean R² = {mean_r2:.3f}")

def main():
    """Main function to run the comprehensive analysis."""
    # Initialize analysis
    analysis = ComprehensiveRegulationAnalysis()
    
    # Run complete analysis
    analysis.run_complete_analysis()

if __name__ == "__main__":
    main()
