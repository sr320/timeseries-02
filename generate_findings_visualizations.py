#!/usr/bin/env python3
"""
Generate comprehensive visualizations for the Context-Dependent Findings Report
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def load_data():
    """Load the analysis results"""
    base_path = Path("code/output/context_dependent_analysis")
    
    # Load main results
    methylation_mirna = pd.read_csv(base_path / "methylation_mirna_context.csv")
    lncrna_mirna = pd.read_csv(base_path / "lncrna_mirna_context.csv")
    multi_way = pd.read_csv(base_path / "multi_way_interactions.csv")
    
    return methylation_mirna, lncrna_mirna, multi_way

def create_overview_chart(methylation_mirna, lncrna_mirna, multi_way):
    """Create overview chart of all findings"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('🧬 Context-Dependent Regulation Analysis - Overview', fontsize=20, fontweight='bold')
    
    # 1. Context-dependent interactions
    context_data = {
        'Methylation-miRNA': [len(methylation_mirna), 
                             len(methylation_mirna[methylation_mirna['context_dependent'] == True])],
        'lncRNA-miRNA': [len(lncrna_mirna), 
                         len(lncrna_mirna[lncrna_mirna['context_dependent'] == True])]
    }
    
    x = list(context_data.keys())
    total = [context_data[k][0] for k in x]
    context_dep = [context_data[k][1] for k in x]
    non_context = [total[i] - context_dep[i] for i in range(len(x))]
    
    ax1.bar(x, total, label='Total Interactions', alpha=0.7, color='lightblue')
    ax1.bar(x, context_dep, label='Context-Dependent', alpha=0.9, color='coral')
    ax1.set_title('Context-Dependent vs Total Interactions', fontweight='bold')
    ax1.set_ylabel('Number of Interactions')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Improvement distributions
    ax2.hist(methylation_mirna['improvement_from_interaction'], bins=30, alpha=0.7, 
             label='Methylation-miRNA', color='skyblue')
    ax2.hist(lncrna_mirna['improvement_from_interaction'], bins=30, alpha=0.7, 
             label='lncRNA-miRNA', color='lightcoral')
    ax2.set_title('Distribution of Interaction Improvements', fontweight='bold')
    ax2.set_xlabel('Improvement from Interaction (R²)')
    ax2.set_ylabel('Frequency')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Multi-way regulation complexity
    regulator_counts = multi_way['n_regulators'].value_counts().sort_index()
    ax3.bar(regulator_counts.index, regulator_counts.values, color='mediumseagreen', alpha=0.8)
    ax3.set_title('Distribution of Regulators per Gene', fontweight='bold')
    ax3.set_xlabel('Number of Regulators')
    ax3.set_ylabel('Number of Genes')
    ax3.grid(True, alpha=0.3)
    
    # 4. Performance improvement
    improvements = multi_way['improvement_from_regulators']
    ax4.hist(improvements, bins=20, color='gold', alpha=0.8, edgecolor='black')
    ax4.axvline(improvements.mean(), color='red', linestyle='--', 
                label=f'Mean: {improvements.mean():.3f}')
    ax4.set_title('Multi-Regulator Model Improvements', fontweight='bold')
    ax4.set_xlabel('Improvement from Regulators (R²)')
    ax4.set_ylabel('Frequency')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Context-Dependent-Findings-Overview.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_context_analysis_charts(methylation_mirna, lncrna_mirna):
    """Create detailed context analysis charts"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('🔍 Detailed Context-Dependent Analysis', fontsize=20, fontweight='bold')
    
    # 1. Context strength comparison
    meth_context = methylation_mirna[methylation_mirna['context_dependent'] == True]['context_strength']
    lnc_context = lncrna_mirna[lncrna_mirna['context_dependent'] == True]['context_strength']
    
    ax1.boxplot([meth_context, lnc_context], labels=['Methylation-miRNA', 'lncRNA-miRNA'])
    ax1.set_title('Context Strength Distribution', fontweight='bold')
    ax1.set_ylabel('Context Strength')
    ax1.grid(True, alpha=0.3)
    
    # 2. Context direction analysis
    meth_direction = methylation_mirna[methylation_mirna['context_dependent'] == True]['context_direction'].value_counts()
    lnc_direction = lncrna_mirna[lncrna_mirna['context_dependent'] == True]['context_direction'].value_counts()
    
    x = np.arange(2)
    width = 0.35
    
    ax2.bar(x - width/2, meth_direction.values, width, label='Methylation-miRNA', alpha=0.8)
    ax2.bar(x + width/2, lnc_direction.values, width, label='lncRNA-miRNA', alpha=0.8)
    ax2.set_title('Context Direction Analysis', fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Negative', 'Positive'])
    ax2.set_ylabel('Count')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Improvement vs Context Strength
    context_meth = methylation_mirna[methylation_mirna['context_dependent'] == True]
    context_lnc = lncrna_mirna[lncrna_mirna['context_dependent'] == True]
    
    ax3.scatter(context_meth['context_strength'], context_meth['improvement_from_interaction'], 
                alpha=0.6, color='skyblue', label='Methylation-miRNA')
    ax3.scatter(context_lnc['context_strength'], context_lnc['improvement_from_interaction'], 
                alpha=0.6, color='lightcoral', label='lncRNA-miRNA')
    ax3.set_title('Context Strength vs Improvement', fontweight='bold')
    ax3.set_xlabel('Context Strength')
    ax3.set_ylabel('Improvement from Interaction (R²)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Correlation patterns
    high_corr_meth = context_meth['corr_high_regulator2'].abs()
    low_corr_meth = context_meth['corr_low_regulator2'].abs()
    
    ax4.scatter(high_corr_meth, low_corr_meth, alpha=0.6, color='mediumseagreen')
    ax4.plot([0, 1], [0, 1], 'r--', alpha=0.8, label='Equal correlation')
    ax4.set_title('High vs Low Context Correlation Patterns', fontweight='bold')
    ax4.set_xlabel('|Correlation| in High Context')
    ax4.set_ylabel('|Correlation| in Low Context')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Context-Dependent-Detailed-Analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_multi_way_charts(multi_way):
    """Create multi-way regulation analysis charts"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('🔄 Multi-Way Regulatory Networks Analysis', fontsize=20, fontweight='bold')
    
    # 1. Regulator type distribution
    regulator_types = []
    for types_str in multi_way['regulator_types']:
        types = eval(types_str)
        for t in types:
            if 'mirna' in t:
                regulator_types.append('miRNA')
            elif 'lncrna' in t:
                regulator_types.append('lncRNA')
            elif 'methylation' in t:
                regulator_types.append('Methylation')
    
    type_counts = pd.Series(regulator_types).value_counts()
    ax1.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%', 
             colors=['lightcoral', 'skyblue', 'lightgreen'])
    ax1.set_title('Distribution of Regulator Types', fontweight='bold')
    
    # 2. Improvement distribution
    ax2.hist(multi_way['improvement_from_regulators'], bins=25, color='gold', alpha=0.8, edgecolor='black')
    ax2.axvline(multi_way['improvement_from_regulators'].mean(), color='red', linestyle='--', 
                label=f'Mean: {multi_way["improvement_from_regulators"].mean():.3f}')
    ax2.set_title('Multi-Regulator Model Improvements', fontweight='bold')
    ax2.set_xlabel('Improvement from Regulators (R²)')
    ax2.set_ylabel('Frequency')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Regulator count vs improvement
    ax3.scatter(multi_way['n_regulators'], multi_way['improvement_from_regulators'], 
                alpha=0.7, color='purple')
    ax3.set_title('Regulator Count vs Improvement', fontweight='bold')
    ax3.set_xlabel('Number of Regulators')
    ax3.set_ylabel('Improvement from Regulators (R²)')
    ax3.grid(True, alpha=0.3)
    
    # 4. Base vs Full model performance
    ax4.scatter(multi_way['r2_base_model'], multi_way['r2_full_model'], 
                alpha=0.7, color='orange')
    ax4.plot([0, 1], [0, 1], 'r--', alpha=0.8, label='Equal performance')
    ax4.set_title('Base Model vs Full Model Performance', fontweight='bold')
    ax4.set_xlabel('Base Model R²')
    ax4.set_ylabel('Full Model R²')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Multi-Way-Regulation-Analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_summary_statistics(methylation_mirna, lncrna_mirna, multi_way):
    """Create summary statistics table"""
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.axis('tight')
    ax.axis('off')
    
    # Prepare data for table
    data = [
        ['Methylation-miRNA Interactions', len(methylation_mirna), 
         len(methylation_mirna[methylation_mirna['context_dependent'] == True]),
         f"{len(methylation_mirna[methylation_mirna['context_dependent'] == True])/len(methylation_mirna)*100:.1f}%",
         f"{methylation_mirna['improvement_from_interaction'].mean():.3f}"],
        ['lncRNA-miRNA Interactions', len(lncrna_mirna),
         len(lncrna_mirna[lncrna_mirna['context_dependent'] == True]),
         f"{len(lncrna_mirna[lncrna_mirna['context_dependent'] == True])/len(lncrna_mirna)*100:.1f}%",
         f"{lncrna_mirna['improvement_from_interaction'].mean():.3f}"],
        ['Multi-Way Regulation', len(multi_way), 
         len(multi_way[multi_way['has_significant_interactions'] == True]),
         f"{len(multi_way[multi_way['has_significant_interactions'] == True])/len(multi_way)*100:.1f}%",
         f"{multi_way['improvement_from_regulators'].mean():.3f}"],
        ['Total Analysis', len(methylation_mirna) + len(lncrna_mirna) + len(multi_way),
         len(methylation_mirna[methylation_mirna['context_dependent'] == True]) + 
         len(lncrna_mirna[lncrna_mirna['context_dependent'] == True]) + 
         len(multi_way[multi_way['has_significant_interactions'] == True]),
         f"{(len(methylation_mirna[methylation_mirna['context_dependent'] == True]) + 
              len(lncrna_mirna[lncrna_mirna['context_dependent'] == True]) + 
              len(multi_way[multi_way['has_significant_interactions'] == True])) / 
              (len(methylation_mirna) + len(lncrna_mirna) + len(multi_way)) * 100:.1f}%",
         "N/A"]
    ]
    
    table = ax.table(cellText=data,
                     colLabels=['Analysis Type', 'Total', 'Significant', 'Percentage', 'Mean Improvement'],
                     cellLoc='center',
                     loc='center',
                     colWidths=[0.25, 0.15, 0.15, 0.15, 0.2])
    
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 2)
    
    # Style the table
    for i in range(len(data) + 1):
        for j in range(5):
            if i == 0:  # Header row
                table[(i, j)].set_facecolor('#4CAF50')
                table[(i, j)].set_text_props(weight='bold', color='white')
            elif i == len(data):  # Total row
                table[(i, j)].set_facecolor('#FF9800')
                table[(i, j)].set_text_props(weight='bold')
            else:
                table[(i, j)].set_facecolor('#E3F2FD')
    
    ax.set_title('📊 Summary Statistics of Context-Dependent Analysis', fontsize=16, fontweight='bold', pad=20)
    plt.savefig('Summary-Statistics-Table.png', dpi=300, bbox_inches='tight')
    plt.show()

def main():
    """Main function to generate all visualizations"""
    print("🧬 Loading context-dependent regulation analysis data...")
    
    try:
        methylation_mirna, lncrna_mirna, multi_way = load_data()
        print(f"✅ Loaded {len(methylation_mirna)} methylation-miRNA interactions")
        print(f"✅ Loaded {len(lncrna_mirna)} lncRNA-miRNA interactions")
        print(f"✅ Loaded {len(multi_way)} multi-way interactions")
        
        print("\n📊 Generating overview chart...")
        create_overview_chart(methylation_mirna, lncrna_mirna, multi_way)
        
        print("🔍 Generating detailed context analysis charts...")
        create_context_analysis_charts(methylation_mirna, lncrna_mirna)
        
        print("🔄 Generating multi-way regulation charts...")
        create_multi_way_charts(multi_way)
        
        print("📋 Generating summary statistics table...")
        create_summary_statistics(methylation_mirna, lncrna_mirna, multi_way)
        
        print("\n🎉 All visualizations generated successfully!")
        print("📁 Files saved:")
        print("   - Context-Dependent-Findings-Overview.png")
        print("   - Context-Dependent-Detailed-Analysis.png")
        print("   - Multi-Way-Regulation-Analysis.png")
        print("   - Summary-Statistics-Table.png")
        
    except Exception as e:
        print(f"❌ Error generating visualizations: {e}")
        print("Please ensure the analysis results are available in code/output/context_dependent_analysis/")

if __name__ == "__main__":
    main()
