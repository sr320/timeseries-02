#!/usr/bin/env python3
"""
Context-Dependent Regulation Analysis

This script identifies context-dependent regulatory interactions where:
1. DNA methylation is predictive of gene expression ONLY when specific miRNAs are highly expressed
2. lncRNA effects on genes depend on miRNA levels
3. Multi-way regulatory interactions that create context-specific effects

Methods used:
- Interaction term analysis
- Conditional correlation analysis
- Multi-variable regression with interaction terms
- Context-specific regulatory network inference
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings('ignore')

class ContextDependentRegulationAnalysis:
    def __init__(self, data_dir="data/cleaned_datasets"):
        """Initialize the context-dependent analysis."""
        self.data_dir = data_dir
        self.datasets = {}
        self.results = {}
        self.load_datasets()
        
    def load_datasets(self):
        """Load all cleaned datasets."""
        print("Loading cleaned datasets for context-dependent analysis...")
        
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
        
    def analyze_context_dependent_regulation(self):
        """Main analysis for context-dependent regulation."""
        print("\n" + "="*80)
        print("CONTEXT-DEPENDENT REGULATION ANALYSIS")
        print("="*80)
        
        # 1. Analyze methylation-gene interactions dependent on miRNA levels
        print("\n1. Analyzing methylation-gene interactions dependent on miRNA levels...")
        methylation_mirna_context = self.analyze_methylation_mirna_context()
        
        # 2. Analyze lncRNA-gene interactions dependent on miRNA levels
        print("\n2. Analyzing lncRNA-gene interactions dependent on miRNA levels...")
        lncrna_mirna_context = self.analyze_lncrna_mirna_context()
        
        # 3. Analyze multi-way regulatory interactions
        print("\n3. Analyzing multi-way regulatory interactions...")
        multi_way_interactions = self.analyze_multi_way_interactions()
        
        # 4. Context-specific regulatory network inference
        print("\n4. Inferring context-specific regulatory networks...")
        context_networks = self.infer_context_specific_networks()
        
        # Store results
        self.results['context_dependent'] = {
            'methylation_mirna_context': methylation_mirna_context,
            'lncrna_mirna_context': lncrna_mirna_context,
            'multi_way_interactions': multi_way_interactions,
            'context_networks': context_networks
        }
        
        return self.results['context_dependent']
    
    def analyze_methylation_mirna_context(self):
        """Analyze how methylation-gene relationships depend on miRNA levels."""
        print("  Analyzing methylation-gene interactions in miRNA context...")
        
        results = []
        
        # Sample a subset for computational efficiency
        n_genes = min(200, len(self.datasets['gene']))
        n_cpgs = min(100, len(self.datasets['methylation']))
        n_mirnas = min(25, len(self.datasets['mirna']))
        
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        sampled_cpgs = np.random.choice(self.datasets['methylation'].index, n_cpgs, replace=False)
        sampled_mirnas = np.random.choice(self.datasets['mirna'].index, n_mirnas, replace=False)
        
        for gene in sampled_genes:
            gene_expression = self.datasets['gene'].loc[gene].values
            
            for cpg in sampled_cpgs:
                methylation_levels = self.datasets['methylation'].loc[cpg].values
                
                for mirna in sampled_mirnas:
                    mirna_levels = self.datasets['mirna'].loc[mirna].values
                    
                    # Analyze context-dependent regulation
                    context_result = self.analyze_triple_interaction(
                        gene_expression, methylation_levels, mirna_levels,
                        gene, cpg, mirna, 'methylation_mirna_gene'
                    )
                    
                    if context_result:
                        results.append(context_result)
        
        return pd.DataFrame(results)
    
    def analyze_lncrna_mirna_context(self):
        """Analyze how lncRNA-gene relationships depend on miRNA levels."""
        print("  Analyzing lncRNA-gene interactions in miRNA context...")
        
        results = []
        
        # Sample a subset for computational efficiency
        n_genes = min(200, len(self.datasets['gene']))
        n_lncrnas = min(100, len(self.datasets['lncrna']))
        n_mirnas = min(25, len(self.datasets['mirna']))
        
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        sampled_lncrnas = np.random.choice(self.datasets['lncrna'].index, n_lncrnas, replace=False)
        sampled_mirnas = np.random.choice(self.datasets['mirna'].index, n_mirnas, replace=False)
        
        for gene in sampled_genes:
            gene_expression = self.datasets['gene'].loc[gene].values
            
            for lncrna in sampled_lncrnas:
                lncrna_levels = self.datasets['lncrna'].loc[lncrna].values
                
                for mirna in sampled_mirnas:
                    mirna_levels = self.datasets['mirna'].loc[mirna].values
                    
                    # Analyze context-dependent regulation
                    context_result = self.analyze_triple_interaction(
                        gene_expression, lncrna_levels, mirna_levels,
                        gene, lncrna, mirna, 'lncrna_mirna_gene'
                    )
                    
                    if context_result:
                        results.append(context_result)
        
        return pd.DataFrame(results)
    
    def analyze_triple_interaction(self, target, regulator1, regulator2, 
                                 target_name, regulator1_name, regulator2_name, interaction_type):
        """Analyze triple interaction between target and two regulators."""
        
        # Prepare data
        data = pd.DataFrame({
            'target': target,
            'regulator1': regulator1,
            'regulator2': regulator2,
            'interaction': regulator1 * regulator2
        })
        
        # Remove any NaN values
        data = data.dropna()
        
        if len(data) < 10:  # Need sufficient data
            return None
        
        # Standardize variables
        scaler = StandardScaler()
        data_scaled = pd.DataFrame(
            scaler.fit_transform(data),
            columns=data.columns
        )
        
        # Model 1: Only regulator1
        model1 = LinearRegression()
        model1.fit(data_scaled[['regulator1']], data_scaled['target'])
        r2_1 = model1.score(data_scaled[['regulator1']], data_scaled['target'])
        
        # Model 2: regulator1 + regulator2
        model2 = LinearRegression()
        model2.fit(data_scaled[['regulator1', 'regulator2']], data_scaled['target'])
        r2_2 = model2.score(data_scaled[['regulator1', 'regulator2']], data_scaled['target'])
        
        # Model 3: regulator1 + regulator2 + interaction
        model3 = LinearRegression()
        model3.fit(data_scaled[['regulator1', 'regulator2', 'interaction']], data_scaled['target'])
        r2_3 = model3.score(data_scaled[['regulator1', 'regulator2', 'interaction']], data_scaled['target'])
        
        # Calculate improvement from interaction term
        improvement_from_interaction = r2_3 - r2_2
        improvement_from_regulator2 = r2_2 - r2_1
        
        # Determine context dependence
        context_dependent = improvement_from_interaction > 0.1  # Significant improvement
        
        # Calculate conditional correlations
        # High regulator2 condition
        high_reg2_mask = data_scaled['regulator2'] > 0.5
        low_reg2_mask = data_scaled['regulator2'] < -0.5
        
        if high_reg2_mask.sum() > 5 and low_reg2_mask.sum() > 5:
            corr_high_reg2, pval_high = pearsonr(
                data_scaled.loc[high_reg2_mask, 'target'],
                data_scaled.loc[high_reg2_mask, 'regulator1']
            )
            corr_low_reg2, pval_low = pearsonr(
                data_scaled.loc[low_reg2_mask, 'target'],
                data_scaled.loc[low_reg2_mask, 'regulator1']
            )
            
            context_strength = abs(corr_high_reg2 - corr_low_reg2)
        else:
            corr_high_reg2 = corr_low_reg2 = context_strength = np.nan
        
        return {
            'interaction_type': interaction_type,
            'target': target_name,
            'regulator1': regulator1_name,
            'regulator2': regulator2_name,
            'r2_regulator1_only': r2_1,
            'r2_regulator1_regulator2': r2_2,
            'r2_with_interaction': r2_3,
            'improvement_from_regulator2': improvement_from_regulator2,
            'improvement_from_interaction': improvement_from_interaction,
            'context_dependent': context_dependent,
            'corr_high_regulator2': corr_high_reg2,
            'corr_low_regulator2': corr_low_reg2,
            'context_strength': context_strength,
            'context_direction': 'positive' if corr_high_reg2 > corr_low_reg2 else 'negative'
        }
    
    def analyze_multi_way_interactions(self):
        """Analyze complex multi-way regulatory interactions."""
        print("  Analyzing multi-way regulatory interactions...")
        
        results = []
        
        # Sample a subset for computational efficiency
        n_genes = min(100, len(self.datasets['gene']))
        n_regulators = min(50, min(len(self.datasets['mirna']), len(self.datasets['lncrna']), len(self.datasets['methylation'])))
        
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        
        for gene in sampled_genes:
            gene_expression = self.datasets['gene'].loc[gene].values
            
            # Get top regulators for this gene
            regulators = self.get_top_regulators(gene, n_regulators)
            
            if len(regulators) >= 3:
                # Analyze multi-way interactions
                multi_way_result = self.analyze_multi_regulator_interaction(
                    gene_expression, regulators, gene
                )
                
                if multi_way_result:
                    results.append(multi_way_result)
        
        return pd.DataFrame(results)
    
    def get_top_regulators(self, gene, n_regulators):
        """Get top regulators for a specific gene."""
        regulators = {}
        
        # Get top miRNA regulators
        for mirna in self.datasets['mirna'].index:
            corr, pval = pearsonr(
                self.datasets['gene'].loc[gene], 
                self.datasets['mirna'].loc[mirna]
            )
            if pval < 0.1:
                regulators[f"mirna_{mirna}"] = self.datasets['mirna'].loc[mirna].values
        
        # Get top lncRNA regulators
        for lncrna in self.datasets['lncrna'].index:
            corr, pval = pearsonr(
                self.datasets['gene'].loc[gene], 
                self.datasets['lncrna'].loc[lncrna]
            )
            if pval < 0.1:
                regulators[f"lncrna_{lncrna}"] = self.datasets['lncrna'].loc[lncrna].values
        
        # Get top methylation regulators
        for cpg in self.datasets['methylation'].index:
            corr, pval = pearsonr(
                self.datasets['gene'].loc[gene], 
                self.datasets['methylation'].loc[cpg]
            )
            if pval < 0.1:
                regulators[f"methylation_{cpg}"] = self.datasets['methylation'].loc[cpg].values
        
        # Sort by correlation strength and take top n
        regulator_corrs = []
        for name, values in regulators.items():
            corr, pval = pearsonr(
                self.datasets['gene'].loc[gene], 
                values
            )
            regulator_corrs.append((name, values, abs(corr)))
        
        regulator_corrs.sort(key=lambda x: x[2], reverse=True)
        return {name: values for name, values, _ in regulator_corrs[:n_regulators]}
    
    def analyze_multi_regulator_interaction(self, target, regulators, target_name):
        """Analyze interactions between multiple regulators."""
        if len(regulators) < 3:
            return None
        
        # Prepare data
        data = pd.DataFrame(regulators)
        data['target'] = target
        
        # Remove any NaN values
        data = data.dropna()
        
        if len(data) < 10:
            return None
        
        # Standardize variables
        scaler = StandardScaler()
        data_scaled = pd.DataFrame(
            scaler.fit_transform(data),
            columns=data.columns
        )
        
        # Get regulator columns
        regulator_cols = [col for col in data_scaled.columns if col != 'target']
        
        # Model 1: Individual regulators only
        model1 = LinearRegression()
        model1.fit(data_scaled[regulator_cols], data_scaled['target'])
        r2_individual = model1.score(data_scaled[regulator_cols], data_scaled['target'])
        
        # Model 2: Add pairwise interactions
        interaction_cols = []
        for i in range(len(regulator_cols)):
            for j in range(i+1, len(regulator_cols)):
                interaction_name = f"interaction_{regulator_cols[i]}_{regulator_cols[j]}"
                data_scaled[interaction_name] = data_scaled[regulator_cols[i]] * data_scaled[regulator_cols[j]]
                interaction_cols.append(interaction_name)
        
        if interaction_cols:
            model2 = LinearRegression()
            model2.fit(data_scaled[regulator_cols + interaction_cols], data_scaled['target'])
            r2_pairwise = model2.score(data_scaled[regulator_cols + interaction_cols], data_scaled['target'])
        else:
            r2_pairwise = r2_individual
        
        # Calculate improvement
        improvement_from_interactions = r2_pairwise - r2_individual
        
        return {
            'target': target_name,
            'n_regulators': len(regulator_cols),
            'r2_individual_regulators': r2_individual,
            'r2_with_pairwise_interactions': r2_pairwise,
            'improvement_from_interactions': improvement_from_interactions,
            'has_significant_interactions': improvement_from_interactions > 0.1
        }
    
    def infer_context_specific_networks(self):
        """Infer regulatory networks that are context-specific."""
        print("  Inferring context-specific regulatory networks...")
        
        # Get context-dependent results
        if 'context_dependent' not in self.results:
            print("    No context-dependent results available. Run analysis first.")
            return {}
        
        context_networks = {}
        
        # 1. High miRNA context network
        high_mirna_context = self.infer_high_mirna_context_network()
        context_networks['high_mirna_context'] = high_mirna_context
        
        # 2. Low miRNA context network
        low_mirna_context = self.infer_low_mirna_context_network()
        context_networks['low_mirna_context'] = low_mirna_context
        
        # 3. High methylation context network
        high_meth_context = self.infer_high_methylation_context_network()
        context_networks['high_methylation_context'] = high_meth_context
        
        return context_networks
    
    def infer_high_mirna_context_network(self):
        """Infer regulatory network when miRNAs are highly expressed."""
        print("    Inferring high miRNA context network...")
        
        # Find samples with high miRNA expression
        mirna_means = self.datasets['mirna'].mean(axis=0)
        high_mirna_samples = mirna_means[mirna_means > mirna_means.quantile(0.75)].index
        
        if len(high_mirna_samples) < 5:
            return {}
        
        # Analyze correlations in high miRNA context
        network = self.analyze_context_specific_correlations(high_mirna_samples, 'high_mirna')
        return network
    
    def infer_low_mirna_context_network(self):
        """Infer regulatory network when miRNAs are lowly expressed."""
        print("    Inferring low miRNA context network...")
        
        # Find samples with low miRNA expression
        mirna_means = self.datasets['mirna'].mean(axis=0)
        low_mirna_samples = mirna_means[mirna_means < mirna_means.quantile(0.25)].index
        
        if len(low_mirna_samples) < 5:
            return {}
        
        # Analyze correlations in low miRNA context
        network = self.analyze_context_specific_correlations(low_mirna_samples, 'low_mirna')
        return network
    
    def infer_high_methylation_context_network(self):
        """Infer regulatory network when methylation is high."""
        print("    Inferring high methylation context network...")
        
        # Find samples with high methylation
        meth_means = self.datasets['methylation'].mean(axis=0)
        high_meth_samples = meth_means[meth_means > meth_means.quantile(0.75)].index
        
        if len(high_meth_samples) < 5:
            return {}
        
        # Analyze correlations in high methylation context
        network = self.analyze_context_specific_correlations(high_meth_samples, 'high_methylation')
        return network
    
    def analyze_context_specific_correlations(self, context_samples, context_name):
        """Analyze correlations in a specific context."""
        network = {
            'context': context_name,
            'n_samples': len(context_samples),
            'gene_mirna_correlations': [],
            'gene_lncrna_correlations': [],
            'gene_methylation_correlations': []
        }
        
        # Sample a subset for efficiency
        n_genes = min(100, len(self.datasets['gene']))
        sampled_genes = np.random.choice(self.datasets['gene'].index, n_genes, replace=False)
        
        # Analyze gene-miRNA correlations in context
        for gene in sampled_genes:
            for mirna in self.datasets['mirna'].index:
                corr, pval = pearsonr(
                    self.datasets['gene'].loc[gene, context_samples],
                    self.datasets['mirna'].loc[mirna, context_samples]
                )
                if pval < 0.05:
                    network['gene_mirna_correlations'].append({
                        'gene': gene,
                        'mirna': mirna,
                        'correlation': corr,
                        'p_value': pval
                    })
        
        # Analyze gene-lncRNA correlations in context
        for gene in sampled_genes:
            for lncrna in self.datasets['lncrna'].index:
                corr, pval = pearsonr(
                    self.datasets['gene'].loc[gene, context_samples],
                    self.datasets['lncrna'].loc[lncrna, context_samples]
                )
                if pval < 0.05:
                    network['gene_lncrna_correlations'].append({
                        'gene': gene,
                        'lncrna': lncrna,
                        'correlation': corr,
                        'p_value': pval
                    })
        
        # Analyze gene-methylation correlations in context
        for gene in sampled_genes:
            for cpg in self.datasets['methylation'].index:
                corr, pval = pearsonr(
                    self.datasets['gene'].loc[gene, context_samples],
                    self.datasets['methylation'].loc[cpg, context_samples]
                )
                if pval < 0.05:
                    network['gene_methylation_correlations'].append({
                        'gene': gene,
                        'cpg': cpg,
                        'correlation': corr,
                        'p_value': pval
                    })
        
        return network
    
    def generate_context_visualizations(self):
        """Generate visualizations for context-dependent analysis."""
        print("\n" + "="*60)
        print("GENERATING CONTEXT-DEPENDENT VISUALIZATIONS")
        print("="*60)
        
        # Create output directory
        import os
        os.makedirs("output/context_dependent_analysis", exist_ok=True)
        
        # 1. Context-dependent correlation plots
        self.plot_context_dependent_correlations()
        
        # 2. Interaction term importance plots
        self.plot_interaction_importance()
        
        # 3. Context-specific network comparisons
        self.plot_context_network_comparisons()
        
        print("✓ All context-dependent visualizations saved to output/context_dependent_analysis/")
    
    def plot_context_dependent_correlations(self):
        """Plot context-dependent correlations."""
        if 'context_dependent' not in self.results:
            print("No context-dependent results to plot.")
            return
        
        # Get methylation-miRNA context results
        meth_mirna = self.results['context_dependent']['methylation_mirna_context']
        
        if not meth_mirna.empty:
            # Plot improvement from interaction terms
            plt.figure(figsize=(12, 8))
            
            plt.subplot(2, 2, 1)
            plt.hist(meth_mirna['improvement_from_interaction'], bins=20, alpha=0.7)
            plt.title('Improvement from Methylation-miRNA Interaction')
            plt.xlabel('R² Improvement')
            plt.ylabel('Frequency')
            
            plt.subplot(2, 2, 2)
            context_strength = meth_mirna['context_strength'].dropna()
            if len(context_strength) > 0:
                plt.hist(context_strength, bins=20, alpha=0.7)
                plt.title('Context Strength (High vs Low miRNA)')
                plt.xlabel('Correlation Difference')
                plt.ylabel('Frequency')
            
            plt.subplot(2, 2, 3)
            high_corr = meth_mirna['corr_high_regulator2'].dropna()
            low_corr = meth_mirna['corr_low_regulator2'].dropna()
            if len(high_corr) > 0 and len(low_corr) > 0:
                plt.scatter(high_corr, low_corr, alpha=0.6)
                plt.plot([-1, 1], [-1, 1], 'r--', alpha=0.5)
                plt.title('Correlation: High vs Low miRNA Context')
                plt.xlabel('Correlation in High miRNA Context')
                plt.ylabel('Correlation in Low miRNA Context')
            
            plt.subplot(2, 2, 4)
            context_dependent = meth_mirna['context_dependent'].value_counts()
            plt.pie(context_dependent.values, labels=context_dependent.index, autopct='%1.1f%%')
            plt.title('Context-Dependent vs Independent')
            
            plt.tight_layout()
            plt.savefig("output/context_dependent_analysis/context_dependent_correlations.png", 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def plot_interaction_importance(self):
        """Plot interaction term importance."""
        if 'context_dependent' not in self.results:
            print("No context-dependent results to plot.")
            return
        
        # Get multi-way interaction results
        multi_way = self.results['context_dependent']['multi_way_interactions']
        
        if not multi_way.empty:
            plt.figure(figsize=(10, 6))
            
            plt.subplot(1, 2, 1)
            plt.hist(multi_way['improvement_from_interactions'], bins=20, alpha=0.7)
            plt.title('Improvement from Multi-Way Interactions')
            plt.xlabel('R² Improvement')
            plt.ylabel('Frequency')
            
            plt.subplot(1, 2, 2)
            significant = multi_way['has_significant_interactions'].value_counts()
            plt.pie(significant.values, labels=significant.index, autopct='%1.1f%%')
            plt.title('Significant Multi-Way Interactions')
            
            plt.tight_layout()
            plt.savefig("output/context_dependent_analysis/interaction_importance.png", 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def plot_context_network_comparisons(self):
        """Plot context-specific network comparisons."""
        if 'context_dependent' not in self.results:
            print("No context-dependent results to plot.")
            return
        
        context_networks = self.results['context_dependent']['context_networks']
        
        if context_networks:
            # Compare network sizes across contexts
            contexts = list(context_networks.keys())
            network_sizes = []
            
            for context in contexts:
                network = context_networks[context]
                total_relationships = (
                    len(network.get('gene_mirna_correlations', [])) +
                    len(network.get('gene_lncrna_correlations', [])) +
                    len(network.get('gene_methylation_correlations', []))
                )
                network_sizes.append(total_relationships)
            
            plt.figure(figsize=(10, 6))
            bars = plt.bar(contexts, network_sizes, alpha=0.7)
            plt.title('Regulatory Network Size by Context')
            plt.ylabel('Number of Regulatory Relationships')
            plt.xticks(rotation=45)
            
            # Add value labels on bars
            for bar, size in zip(bars, network_sizes):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01*max(network_sizes),
                        str(size), ha='center', va='bottom')
            
            plt.tight_layout()
            plt.savefig("output/context_dependent_analysis/context_network_comparisons.png", 
                       dpi=300, bbox_inches='tight')
            plt.close()
    
    def save_context_results(self):
        """Save context-dependent analysis results."""
        print("\n" + "="*60)
        print("SAVING CONTEXT-DEPENDENT RESULTS")
        print("="*60)
        
        import os
        os.makedirs("output/context_dependent_analysis", exist_ok=True)
        
        # Save context-dependent results
        if 'context_dependent' in self.results:
            for analysis_type, results in self.results['context_dependent'].items():
                if isinstance(results, pd.DataFrame) and not results.empty:
                    results.to_csv(f"output/context_dependent_analysis/{analysis_type}.csv")
                    print(f"Saved {analysis_type}: {len(results)} results")
                elif isinstance(results, dict):
                    # Save context networks
                    for context_name, context_data in results.items():
                        if isinstance(context_data, dict):
                            # Save each correlation type
                            for corr_type, corr_data in context_data.items():
                                if isinstance(corr_data, list) and corr_data:
                                    corr_df = pd.DataFrame(corr_data)
                                    corr_df.to_csv(f"output/context_dependent_analysis/{context_name}_{corr_type}.csv")
                                    print(f"Saved {context_name}_{corr_type}: {len(corr_data)} correlations")
        
        print("✓ All context-dependent results saved to output/context_dependent_analysis/")
    
    def print_context_summary(self):
        """Print summary of context-dependent analysis."""
        print("\n" + "="*60)
        print("CONTEXT-DEPENDENT ANALYSIS SUMMARY")
        print("="*60)
        
        if 'context_dependent' not in self.results:
            print("No context-dependent results available.")
            return
        
        # Summary of methylation-miRNA context
        meth_mirna = self.results['context_dependent']['methylation_mirna_context']
        if not meth_mirna.empty:
            print(f"\nMETHYLATION-MIRNA CONTEXT ANALYSIS:")
            print(f"  Total interactions analyzed: {len(meth_mirna)}")
            print(f"  Context-dependent interactions: {meth_mirna['context_dependent'].sum()}")
            print(f"  Mean improvement from interaction: {meth_mirna['improvement_from_interaction'].mean():.3f}")
            print(f"  Mean context strength: {meth_mirna['context_strength'].mean():.3f}")
        
        # Summary of multi-way interactions
        multi_way = self.results['context_dependent']['multi_way_interactions']
        if not multi_way.empty:
            print(f"\nMULTI-WAY INTERACTION ANALYSIS:")
            print(f"  Total genes analyzed: {len(multi_way)}")
            print(f"  Genes with significant interactions: {multi_way['has_significant_interactions'].sum()}")
            print(f"  Mean improvement from interactions: {multi_way['improvement_from_interactions'].mean():.3f}")
        
        # Summary of context networks
        context_networks = self.results['context_dependent']['context_networks']
        if context_networks:
            print(f"\nCONTEXT-SPECIFIC NETWORKS:")
            for context_name, network in context_networks.items():
                total_relationships = (
                    len(network.get('gene_mirna_correlations', [])) +
                    len(network.get('gene_lncrna_correlations', [])) +
                    len(network.get('gene_methylation_correlations', []))
                )
                print(f"  {context_name}: {total_relationships} regulatory relationships")
    
    def run_complete_context_analysis(self):
        """Run the complete context-dependent analysis."""
        print("="*80)
        print("CONTEXT-DEPENDENT REGULATION ANALYSIS")
        print("="*80)
        print("This analysis will identify context-dependent regulatory interactions")
        print("where the effect of one regulator depends on the level of another.")
        print("="*80)
        
        # Run context-dependent analysis
        self.analyze_context_dependent_regulation()
        
        # Generate visualizations
        self.generate_context_visualizations()
        
        # Save results
        self.save_context_results()
        
        # Print summary
        self.print_context_summary()
        
        print("\n" + "="*80)
        print("CONTEXT-DEPENDENT ANALYSIS COMPLETED SUCCESSFULLY!")
        print("="*80)

def main():
    """Main function to run the context-dependent analysis."""
    # Initialize analysis
    analysis = ContextDependentRegulationAnalysis()
    
    # Run complete context-dependent analysis
    analysis.run_complete_context_analysis()

if __name__ == "__main__":
    main()
