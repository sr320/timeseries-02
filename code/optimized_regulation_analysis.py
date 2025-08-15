#!/usr/bin/env python3
"""
OPTIMIZED Comprehensive Analysis: How miRNA, lncRNA, and DNA methylation influence gene expression

This optimized script utilizes:
- Parallel processing across 48 CPU cores
- Memory-efficient batch operations using 247GB RAM
- Vectorized numpy/pandas operations
- Chunked processing for large datasets
- Concurrent data loading and processing
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
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import partial
import os
import gc
import time
from typing import Dict, List, Tuple, Any
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class OptimizedRegulationAnalysis:
    def __init__(self, data_dir="data/cleaned_datasets", n_jobs=None):
        """Initialize the optimized analysis with parallel processing capabilities."""
        self.data_dir = data_dir
        self.datasets = {}
        self.results = {}
        
        # Set number of jobs for parallel processing
        if n_jobs is None:
            self.n_jobs = min(48, mp.cpu_count())  # Use all available cores
        else:
            self.n_jobs = n_jobs
            
        print(f"🚀 Initializing optimized analysis with {self.n_jobs} parallel workers")
        print(f"💾 Available RAM: {self._get_available_ram():.1f} GB")
        
        self.load_datasets()
        
    def _get_available_ram(self):
        """Get available RAM in GB."""
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemAvailable:'):
                        return int(line.split()[1]) / (1024**2)
        except:
            pass
        return 200  # Default fallback
        
    def load_datasets(self):
        """Load all cleaned datasets with memory optimization."""
        print("📂 Loading cleaned datasets with memory optimization...")
        
        # Load datasets in parallel
        with ThreadPoolExecutor(max_workers=4) as executor:
            future_gene = executor.submit(pd.read_csv, f"{self.data_dir}/gene_counts_cleaned.csv", index_col=0)
            future_lncrna = executor.submit(pd.read_csv, f"{self.data_dir}/lncrna_counts_cleaned.csv", index_col=0)
            future_mirna = executor.submit(pd.read_csv, f"{self.data_dir}/mirna_counts_cleaned.csv", index_col=0)
            future_methylation = executor.submit(pd.read_csv, f"{self.data_dir}/wgbs_counts_cleaned.csv", index_col=0)
            
            self.datasets['gene'] = future_gene.result()
            self.datasets['lncrna'] = future_lncrna.result()
            self.datasets['mirna'] = future_mirna.result()
            self.datasets['methylation'] = future_methylation.result()
        
        print(f"✅ Loaded gene expression: {self.datasets['gene'].shape}")
        print(f"✅ Loaded lncRNA expression: {self.datasets['lncrna'].shape}")
        print(f"✅ Loaded miRNA expression: {self.datasets['mirna'].shape}")
        print(f"✅ Loaded DNA methylation: {self.datasets['methylation'].shape}")
        
        # Verify sample alignment
        self.verify_sample_alignment()
        
        # Pre-compute correlation matrices for memory efficiency
        self._precompute_correlations()
        
    def _precompute_correlations(self):
        """Pre-compute correlation matrices to avoid repeated calculations."""
        print("🔢 Pre-computing correlation matrices...")
        
        # Convert to numpy arrays for faster operations
        self.gene_array = self.datasets['gene'].values
        self.lncrna_array = self.datasets['lncrna'].values
        self.mirna_array = self.datasets['mirna'].values
        self.methylation_array = self.datasets['methylation'].values
        
        # Pre-compute sample correlations (vectorized)
        self.sample_correlations = {}
        for data_type, data_array in [('gene', self.gene_array), ('lncrna', self.lncrna_array), 
                                     ('mirna', self.mirna_array), ('methylation', self.methylation_array)]:
            # Vectorized correlation calculation
            corr_matrix = np.corrcoef(data_array.T)
            self.sample_correlations[data_type] = corr_matrix
            
        print("✅ Correlation matrices pre-computed")
        
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
        
    def vectorized_correlation_analysis(self):
        """Perform vectorized correlation analysis between data types."""
        print("\n" + "="*60)
        print("VECTORIZED CORRELATION ANALYSIS")
        print("="*60)
        
        # Use vectorized operations for cross-correlations
        results = {}
        
        # Gene vs lncRNA correlations (vectorized)
        print("🔬 Analyzing gene-lncRNA correlations (vectorized)...")
        gene_lncrna_corr = self._vectorized_cross_correlations(
            self.gene_array, self.lncrna_array, 'gene_lncrna'
        )
        results['gene_lncrna'] = gene_lncrna_corr
        
        # Gene vs miRNA correlations (vectorized)
        print("🔬 Analyzing gene-miRNA correlations (vectorized)...")
        gene_mirna_corr = self._vectorized_cross_correlations(
            self.gene_array, self.mirna_array, 'gene_mirna'
        )
        results['gene_mirna'] = gene_mirna_corr
        
        # Gene vs methylation correlations (vectorized)
        print("🔬 Analyzing gene-methylation correlations (vectorized)...")
        gene_meth_corr = self._vectorized_cross_correlations(
            self.gene_array, self.methylation_array, 'gene_methylation'
        )
        results['gene_methylation'] = gene_meth_corr
        
        # lncRNA vs miRNA correlations
        print("🔬 Analyzing lncRNA-miRNA correlations (vectorized)...")
        lncrna_mirna_corr = self._vectorized_cross_correlations(
            self.lncrna_array, self.mirna_array, 'lncrna_mirna'
        )
        results['lncrna_mirna'] = lncrna_mirna_corr
        
        self.results['correlations'] = results
        print("✅ Vectorized correlation analysis completed")
        
    def _vectorized_cross_correlations(self, data1: np.ndarray, data2: np.ndarray, 
                                     analysis_type: str) -> Dict[str, Any]:
        """Calculate cross-correlations using vectorized operations."""
        n_features1, n_features2 = data1.shape[0], data2.shape[0]
        
        # Use chunked processing for very large matrices to manage memory
        chunk_size = 1000  # Process 1000 features at a time
        all_correlations = []
        
        for i in range(0, n_features1, chunk_size):
            end_i = min(i + chunk_size, n_features1)
            chunk_correlations = []
            
            for j in range(0, n_features2, chunk_size):
                end_j = min(j + chunk_size, n_features2)
                
                # Vectorized correlation calculation for this chunk
                chunk_data1 = data1[i:end_i]
                chunk_data2 = data2[j:end_j]
                
                # Calculate correlations for all pairs in this chunk
                correlations = np.corrcoef(chunk_data1, chunk_data2)[:chunk_data1.shape[0], chunk_data1.shape[0]:]
                
                # Store significant correlations
                for idx1 in range(correlations.shape[0]):
                    for idx2 in range(correlations.shape[1]):
                        corr_val = correlations[idx1, idx2]
                        if not np.isnan(corr_val) and abs(corr_val) > 0.3:  # Filter significant correlations
                            chunk_correlations.append({
                                'feature1_idx': i + idx1,
                                'feature2_idx': j + idx2,
                                'correlation': corr_val
                            })
            
            all_correlations.extend(chunk_correlations)
            
            # Memory cleanup
            gc.collect()
        
        # Convert to DataFrame for easier analysis
        if all_correlations:
            corr_df = pd.DataFrame(all_correlations)
            corr_df['feature1'] = corr_df['feature1_idx'].apply(lambda x: self.datasets[analysis_type.split('_')[0]].index[x])
            corr_df['feature2'] = corr_df['feature2_idx'].apply(lambda x: self.datasets[analysis_type.split('_')[1]].index[x])
            corr_df = corr_df.drop(['feature1_idx', 'feature2_idx'], axis=1)
        else:
            corr_df = pd.DataFrame()
        
        return {
            'correlations': corr_df,
            'n_significant': len(all_correlations),
            'mean_correlation': np.mean([c['correlation'] for c in all_correlations]) if all_correlations else 0
        }
        
    def parallel_regulatory_network_analysis(self):
        """Perform regulatory network analysis using parallel processing."""
        print("\n" + "="*60)
        print("PARALLEL REGULATORY NETWORK ANALYSIS")
        print("="*60)
        
        # Split genes into chunks for parallel processing
        gene_chunks = np.array_split(self.datasets['gene'].index, self.n_jobs)
        
        print(f"🔄 Processing {len(self.datasets['gene'].index)} genes in {len(gene_chunks)} parallel chunks...")
        
        # Process chunks in parallel
        with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            # Submit all chunks for processing
            future_to_chunk = {
                executor.submit(self._process_gene_chunk, chunk): i 
                for i, chunk in enumerate(gene_chunks)
            }
            
            # Collect results as they complete
            all_networks = {'mirna_gene': [], 'lncrna_gene': [], 'methylation_gene': []}
            
            for future in as_completed(future_to_chunk):
                chunk_idx = future_to_chunk[future]
                try:
                    chunk_results = future.result()
                    for network_type in all_networks:
                        all_networks[network_type].extend(chunk_results[network_type])
                    print(f"✅ Chunk {chunk_idx + 1}/{len(gene_chunks)} completed")
                except Exception as exc:
                    print(f"❌ Chunk {chunk_idx + 1} generated an exception: {exc}")
        
        # Combine results
        self.results['regulatory_networks'] = {}
        for network_type, networks in all_networks.items():
            if networks:
                self.results['regulatory_networks'][network_type] = pd.DataFrame(networks)
                print(f"✅ {network_type}: {len(networks)} regulatory relationships identified")
            else:
                self.results['regulatory_networks'][network_type] = pd.DataFrame()
                
        print("✅ Parallel regulatory network analysis completed")
        
    def _process_gene_chunk(self, gene_chunk: List[str]) -> Dict[str, List]:
        """Process a chunk of genes for regulatory network analysis."""
        networks = {'mirna_gene': [], 'lncrna_gene': [], 'methylation_gene': []}
        
        for gene in gene_chunk:
            gene_expression = self.datasets['gene'].loc[gene].values
            
            # miRNA regulators (top correlations)
            mirna_corrs = []
            for mirna in self.datasets['mirna'].index:
                corr, pval = pearsonr(gene_expression, self.datasets['mirna'].loc[mirna])
                if pval < 0.1:
                    mirna_corrs.append((mirna, corr, pval))
            
            # Sort and take top correlations
            mirna_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
            for mirna, corr, pval in mirna_corrs[:5]:
                networks['mirna_gene'].append({
                    'gene': gene, 'mirna': mirna, 'correlation': corr, 'p_value': pval
                })
            
            # lncRNA regulators
            lncrna_corrs = []
            for lncrna in self.datasets['lncrna'].index:
                corr, pval = pearsonr(gene_expression, self.datasets['lncrna'].loc[lncrna])
                if pval < 0.1:
                    lncrna_corrs.append((lncrna, corr, pval))
            
            lncrna_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
            for lncrna, corr, pval in lncrna_corrs[:5]:
                networks['lncrna_gene'].append({
                    'gene': gene, 'lncrna': lncrna, 'correlation': corr, 'p_value': pval
                })
            
            # Methylation regulators
            meth_corrs = []
            for cpg in self.datasets['methylation'].index:
                corr, pval = pearsonr(gene_expression, self.datasets['methylation'].loc[cpg])
                if pval < 0.1:
                    meth_corrs.append((cpg, corr, pval))
            
            meth_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
            for cpg, corr, pval in meth_corrs[:5]:
                networks['methylation_gene'].append({
                    'gene': gene, 'cpg': cpg, 'correlation': corr, 'p_value': pval
                })
        
        return networks
        
    def parallel_predictive_modeling(self):
        """Train predictive models using parallel processing."""
        print("\n" + "="*60)
        print("PARALLEL PREDICTIVE MODELING")
        print("="*60)
        
        # Prepare integrated data for modeling
        print("🔧 Preparing integrated dataset for modeling...")
        integrated_data = self._prepare_integrated_dataset()
        
        # Split genes into chunks for parallel processing
        gene_chunks = np.array_split(list(integrated_data.keys()), self.n_jobs)
        
        print(f"🔄 Training models for {len(integrated_data)} genes in {len(gene_chunks)} parallel chunks...")
        
        # Train models in parallel
        with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
            future_to_chunk = {
                executor.submit(self._train_models_chunk, chunk, integrated_data): i 
                for i, chunk in enumerate(gene_chunks)
            }
            
            all_models = {'linear': {}, 'ridge': {}, 'lasso': {}, 'random_forest': {}}
            
            for future in as_completed(future_to_chunk):
                chunk_idx = future_to_chunk[future]
                try:
                    chunk_results = future.result()
                    for model_type in all_models:
                        all_models[model_type].update(chunk_results[model_type])
                    print(f"✅ Model training chunk {chunk_idx + 1}/{len(gene_chunks)} completed")
                except Exception as exc:
                    print(f"❌ Model training chunk {chunk_idx + 1} generated an exception: {exc}")
        
        self.results['predictive_models'] = all_models
        print("✅ Parallel predictive modeling completed")
        
    def _prepare_integrated_dataset(self) -> Dict[str, pd.DataFrame]:
        """Prepare integrated dataset for modeling."""
        integrated_data = {}
        
        print("🔗 Integrating regulatory factors for each gene...")
        
        # Process genes in batches for memory efficiency
        batch_size = 1000
        for i in range(0, len(self.datasets['gene'].index), batch_size):
            batch_genes = self.datasets['gene'].index[i:i+batch_size]
            
            for gene in batch_genes:
                gene_expression = self.datasets['gene'].loc[gene].values
                regulators = {}
                
                # Top miRNA regulators
                mirna_corrs = []
                for mirna in self.datasets['mirna'].index:
                    corr, pval = pearsonr(gene_expression, self.datasets['mirna'].loc[mirna])
                    if pval < 0.1:
                        mirna_corrs.append((mirna, corr, pval))
                
                mirna_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
                for mirna, corr, pval in mirna_corrs[:5]:
                    regulators[f"mirna_{mirna}"] = self.datasets['mirna'].loc[mirna].values
                
                # Top lncRNA regulators
                lncrna_corrs = []
                for lncrna in self.datasets['lncrna'].index:
                    corr, pval = pearsonr(gene_expression, self.datasets['lncrna'].loc[lncrna])
                    if pval < 0.1:
                        lncrna_corrs.append((lncrna, corr, pval))
                
                lncrna_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
                for lncrna, corr, pval in lncrna_corrs[:5]:
                    regulators[f"lncrna_{lncrna}"] = self.datasets['lncrna'].loc[lncrna].values
                
                # Top methylation regulators
                meth_corrs = []
                for cpg in self.datasets['methylation'].index:
                    corr, pval = pearsonr(gene_expression, self.datasets['methylation'].loc[cpg])
                    if pval < 0.1:
                        meth_corrs.append((cpg, corr, pval))
                
                meth_corrs.sort(key=lambda x: abs(x[1]), reverse=True)
                for cpg, corr, pval in meth_corrs[:5]:
                    regulators[f"methylation_{cpg}"] = self.datasets['methylation'].loc[cpg].values
                
                # Create feature matrix
                if regulators:
                    feature_matrix = pd.DataFrame(regulators)
                    feature_matrix['gene_expression'] = gene_expression
                    integrated_data[gene] = feature_matrix
            
            # Memory cleanup after each batch
            gc.collect()
        
        return integrated_data
        
    def _train_models_chunk(self, gene_chunk: List[str], integrated_data: Dict[str, pd.DataFrame]) -> Dict[str, Dict]:
        """Train models for a chunk of genes."""
        models = {'linear': {}, 'ridge': {}, 'lasso': {}, 'random_forest': {}}
        
        for gene in gene_chunk:
            if gene in integrated_data:
                data = integrated_data[gene]
                if len(data) > 10:  # Need sufficient data
                    X = data.drop('gene_expression', axis=1)
                    y = data['gene_expression']
                    
                    # Split data
                    X_train, X_test, y_train, y_test = train_test_split(
                        X, y, test_size=0.3, random_state=42
                    )
                    
                    # Train Linear Regression
                    try:
                        linear_model = LinearRegression()
                        linear_model.fit(X_train, y_train)
                        y_pred = linear_model.predict(X_test)
                        r2 = r2_score(y_test, y_pred)
                        mse = mean_squared_error(y_test, y_pred)
                        
                        models['linear'][gene] = {
                            'r2_score': r2,
                            'mse': mse,
                            'feature_importance': dict(zip(X.columns, linear_model.coef_))
                        }
                    except:
                        pass
                    
                    # Train Ridge Regression
                    try:
                        ridge_model = Ridge(alpha=1.0)
                        ridge_model.fit(X_train, y_train)
                        y_pred = ridge_model.predict(X_test)
                        r2 = r2_score(y_test, y_pred)
                        mse = mean_squared_error(y_test, y_pred)
                        
                        models['ridge'][gene] = {
                            'r2_score': r2,
                            'mse': mse,
                            'feature_importance': dict(zip(X.columns, ridge_model.coef_))
                        }
                    except:
                        pass
                    
                    # Train Lasso Regression
                    try:
                        lasso_model = Lasso(alpha=0.1)
                        lasso_model.fit(X_train, y_train)
                        y_pred = lasso_model.predict(X_test)
                        r2 = r2_score(y_test, y_pred)
                        mse = mean_squared_error(y_test, y_pred)
                        
                        models['lasso'][gene] = {
                            'r2_score': r2,
                            'mse': mse,
                            'feature_importance': dict(zip(X.columns, lasso_model.coef_))
                        }
                    except:
                        pass
                    
                    # Train Random Forest
                    try:
                        rf_model = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=1)
                        rf_model.fit(X_train, y_train)
                        y_pred = rf_model.predict(X_test)
                        r2 = r2_score(y_test, y_pred)
                        mse = mean_squared_error(y_test, y_pred)
                        
                        models['random_forest'][gene] = {
                            'r2_score': r2,
                            'mse': mse,
                            'feature_importance': dict(zip(X.columns, rf_model.feature_importances_))
                        }
                    except:
                        pass
        
        return models
        
    def vectorized_time_series_analysis(self):
        """Perform vectorized time series analysis."""
        print("\n" + "="*60)
        print("VECTORIZED TIME SERIES ANALYSIS")
        print("="*60)
        
        results = {}
        
        for condition in self.conditions:
            results[condition] = {}
            condition_samples = [s for s in self.samples if f"ACR-{condition}-" in s]
            
            for data_type in ['gene', 'lncrna', 'mirna', 'methylation']:
                data = self.datasets[data_type]
                condition_data = data[condition_samples]
                
                # Vectorized time point analysis
                time_trends = {}
                for tp in self.time_points:
                    tp_samples = [s for s in condition_samples if s.endswith(tp)]
                    if tp_samples:
                        # Vectorized mean calculation
                        tp_data = condition_data[tp_samples]
                        time_trends[tp] = {
                            'mean_expression': tp_data.mean(axis=1).to_dict(),
                            'std_expression': tp_data.std(axis=1).to_dict(),
                            'n_samples': len(tp_samples)
                        }
                
                results[condition][data_type] = time_trends
        
        self.results['time_series'] = results
        print("✅ Vectorized time series analysis completed")
        
    def generate_visualizations(self):
        """Generate comprehensive visualizations."""
        print("\n" + "="*60)
        print("GENERATING VISUALIZATIONS")
        print("="*60)
        
        # Create output directory
        os.makedirs("output/regulation_analysis", exist_ok=True)
        
        # Generate plots in parallel
        with ThreadPoolExecutor(max_workers=4) as executor:
            future_corr = executor.submit(self.plot_correlation_heatmaps)
            future_time = executor.submit(self.plot_time_series)
            future_networks = executor.submit(self.plot_regulatory_networks)
            future_models = executor.submit(self.plot_model_performance)
            
            # Wait for all visualizations to complete
            future_corr.result()
            future_time.result()
            future_networks.result()
            future_models.result()
        
        print("✅ All visualizations saved to output/regulation_analysis/")
        
    def plot_correlation_heatmaps(self):
        """Plot correlation heatmaps using pre-computed correlations."""
        # Use pre-computed sample correlations
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.ravel()
        
        for i, (data_type, corr_matrix) in enumerate(self.sample_correlations.items()):
            # Plot first 20x20 for visualization
            sns.heatmap(corr_matrix[:20, :20], cmap='RdBu_r', center=0, ax=axes[i])
            axes[i].set_title(f'{data_type.title()} Sample Correlations')
        
        plt.tight_layout()
        plt.savefig("output/regulation_analysis/sample_correlation_heatmaps.png", dpi=300, bbox_inches='tight')
        plt.close()
        
    def plot_time_series(self):
        """Plot time series analysis."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.ravel()
        
        for i, data_type in enumerate(['gene', 'lncrna', 'mirna', 'methylation']):
            data = self.datasets[data_type]
            
            # Calculate mean expression per time point per condition
            time_means = {}
            for condition in self.conditions:
                condition_samples = [s for s in self.samples if f"ACR-{condition}-" in s]
                condition_data = data[condition_samples]
                
                for tp in self.time_points:
                    tp_samples = [s for s in condition_samples if s.endswith(tp)]
                    if tp_samples:
                        tp_data = condition_data[tp_samples]
                        time_means[f"{condition}_{tp}"] = tp_data.mean(axis=1).mean()
            
            # Plot time series
            time_points = list(time_means.keys())
            values = list(time_means.values())
            
            axes[i].plot(range(len(time_points)), values, 'o-', linewidth=2, markersize=6)
            axes[i].set_title(f'{data_type.title()} Expression Over Time')
            axes[i].set_xlabel('Time Point')
            axes[i].set_ylabel('Mean Expression')
            axes[i].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig("output/regulation_analysis/time_series_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        
    def plot_regulatory_networks(self):
        """Plot regulatory network visualizations."""
        if 'regulatory_networks' not in self.results:
            return
            
        # Plot network density
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        for i, (network_type, network_data) in enumerate(self.results['regulatory_networks'].items()):
            if isinstance(network_data, pd.DataFrame) and not network_data.empty:
                # Plot correlation distribution
                axes[i].hist(network_data['correlation'], bins=30, alpha=0.7, edgecolor='black')
                axes[i].set_title(f'{network_type.replace("_", " ").title()} Network')
                axes[i].set_xlabel('Correlation Coefficient')
                axes[i].set_ylabel('Frequency')
                axes[i].axvline(x=0, color='red', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig("output/regulation_analysis/regulation_type_distributions.png", dpi=300, bbox_inches='tight')
        plt.close()
        
    def plot_model_performance(self):
        """Plot model performance comparisons."""
        if 'predictive_models' not in self.results:
            return
            
        # Collect performance metrics
        model_performance = {}
        for model_type, gene_results in self.results['predictive_models'].items():
            if gene_results:
                r2_scores = [r['r2_score'] for r in gene_results.values() if 'r2_score' in r]
                if r2_scores:
                    model_performance[model_type] = r2_scores
        
        if model_performance:
            # Create box plot
            plt.figure(figsize=(10, 6))
            plt.boxplot(model_performance.values(), labels=model_performance.keys())
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
        
        # Create output directory
        os.makedirs("output/regulation_analysis", exist_ok=True)
        
        # Save regulatory networks
        if 'regulatory_networks' in self.results:
            for network_type, network_data in self.results['regulatory_networks'].items():
                if isinstance(network_data, pd.DataFrame) and not network_data.empty:
                    network_data.to_csv(f"output/regulation_analysis/{network_type}.csv", index=False)
        
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
                if gene_results:
                    r2_scores = [r['r2_score'] for r in gene_results.values() if 'r2_score' in r]
                    mse_scores = [r['mse'] for r in gene_results.values() if 'mse' in r]
                    
                    model_summary[model_type] = {
                        'mean_r2': np.mean(r2_scores) if r2_scores else 0,
                        'mean_mse': np.mean(mse_scores) if mse_scores else 0,
                        'n_genes': len(gene_results)
                    }
            
            pd.DataFrame(model_summary).T.to_csv("output/regulation_analysis/model_performance_summary.csv")
        
        print("✅ All results saved to output/regulation_analysis/")
        
    def run_complete_analysis(self):
        """Run the complete optimized analysis."""
        start_time = time.time()
        
        print("="*80)
        print("🚀 OPTIMIZED COMPREHENSIVE REGULATION ANALYSIS")
        print("="*80)
        print("This optimized analysis utilizes:")
        print(f"  • {self.n_jobs} parallel CPU cores")
        print(f"  • Vectorized numpy/pandas operations")
        print(f"  • Memory-efficient batch processing")
        print(f"  • Concurrent data loading and processing")
        print("="*80)
        
        # Run all analyses
        self.vectorized_correlation_analysis()
        self.vectorized_time_series_analysis()
        self.parallel_regulatory_network_analysis()
        self.parallel_predictive_modeling()
        
        # Generate visualizations
        self.generate_visualizations()
        
        # Save results
        self.save_results()
        
        # Print summary
        self.print_analysis_summary()
        
        total_time = time.time() - start_time
        print(f"\n⏱️  Total analysis time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
        print("\n" + "="*80)
        print("🎉 OPTIMIZED ANALYSIS COMPLETED SUCCESSFULLY!")
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
                if 'n_significant' in corr_data:
                    print(f"  {corr_type}: {corr_data['n_significant']} significant correlations")
                    print(f"    Mean correlation: {corr_data['mean_correlation']:.3f}")
        
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
                    r2_scores = [r['r2_score'] for r in gene_results.values() if 'r2_score' in r]
                    if r2_scores:
                        mean_r2 = np.mean(r2_scores)
                        print(f"  {model_type}: Mean R² = {mean_r2:.3f} ({len(gene_results)} genes)")

def main():
    """Main function to run the optimized comprehensive analysis."""
    # Initialize optimized analysis
    analysis = OptimizedRegulationAnalysis()
    
    # Run complete analysis
    analysis.run_complete_analysis()

if __name__ == "__main__":
    main()
