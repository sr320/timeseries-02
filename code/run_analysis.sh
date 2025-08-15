#!/bin/bash

# Gene Regulation Analysis Runner Script
# This script sets up the environment and runs the analysis pipeline

echo "=========================================="
echo "Gene Regulation Analysis Pipeline"
echo "=========================================="

# Check if we're in the right directory
if [ ! -d "data" ] || [ ! -d "code" ] || [ ! -d "output" ]; then
    echo "Error: Please run this script from the project root directory"
    echo "Expected structure:"
    echo "  ./data/     (input data files)"
    echo "  ./code/     (analysis scripts)"
    echo "  ./output/   (output directory)"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p output

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    exit 1
fi

# Check if required packages are installed
echo "Checking Python packages..."
python3 -c "
import pandas, numpy, matplotlib, seaborn, scipy, sklearn
print('All required packages are available')
" 2>/dev/null

if [ $? -ne 0 ]; then
    echo "Installing required packages..."
    pip3 install -r code/requirements.txt
    if [ $? -ne 0 ]; then
        echo "Error: Failed to install required packages"
        exit 1
    fi
fi

# Check if data files exist
echo "Checking data files..."
required_files=(
    "data/apul-gene_count_matrix.csv"
    "data/Apul_miRNA_counts_formatted.txt"
    "data/Apul_lncRNA_counts_filtered.txt"
    "data/merged-WGBS-CpG-counts_filtered.csv"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required data file not found: $file"
        exit 1
    fi
done

echo "All data files found"

# Run the simple correlation analysis first (more robust)
echo ""
echo "Running Simple Correlation Analysis..."
echo "====================================="
cd code
python3 simple_correlation_analysis.py
correlation_exit=$?

if [ $correlation_exit -eq 0 ]; then
    echo ""
    echo "Correlation analysis completed successfully!"
else
    echo ""
    echo "Correlation analysis failed with exit code $correlation_exit"
fi

# Try to run the comprehensive analysis if correlation analysis succeeded
if [ $correlation_exit -eq 0 ]; then
    echo ""
    echo "Running Comprehensive Regression Analysis..."
    echo "=========================================="
    python3 analyze_gene_regulation.py
    regression_exit=$?
    
    if [ $regression_exit -eq 0 ]; then
        echo ""
        echo "Comprehensive analysis completed successfully!"
    else
        echo ""
        echo "Comprehensive analysis failed with exit code $regression_exit"
        echo "Correlation analysis results are still available in output/"
    fi
fi

cd ..

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Output files are available in the output/ directory:"
echo ""

# List output files
if [ -d "output" ]; then
    echo "Generated files:"
    ls -la output/
    echo ""
    echo "Summary files:"
    if [ -f "output/correlation_summary.txt" ]; then
        echo "  - correlation_summary.txt (correlation analysis results)"
    fi
    if [ -f "output/analysis_summary.txt" ]; then
        echo "  - analysis_summary.txt (comprehensive analysis results)"
    fi
    echo ""
    echo "Visualization files:"
    ls -la output/*.png 2>/dev/null | head -10
    echo ""
    echo "Data files:"
    ls -la output/*.csv 2>/dev/null | head -10
fi

echo ""
echo "Check the output/ directory for detailed results and visualizations."
echo "=========================================="
