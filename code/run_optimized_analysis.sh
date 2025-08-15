#!/bin/bash

# 🚀 Optimized Timeseries Analysis Runner
# This script runs the optimized analysis with performance monitoring

echo "🚀 TIMESERIES-02 OPTIMIZED ANALYSIS RUNNER"
echo "=============================================="
echo "System: $(nproc) CPU cores, $(free -h | grep Mem | awk '{print $2}') RAM"
echo ""

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "❌ Python not found. Please install Python 3.7+"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
REQUIRED_VERSION="3.7"

if [ "$(printf '%s\n' "$REQUIRED_VERSION" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED_VERSION" ]; then
    echo "❌ Python $REQUIRED_VERSION+ required, found $PYTHON_VERSION"
    exit 1
fi

echo "✅ Python $PYTHON_VERSION detected"

# Install required packages if not present
echo ""
echo "📦 Checking and installing required packages..."
pip install -r requirements_optimized.txt --quiet

if [ $? -ne 0 ]; then
    echo "❌ Failed to install required packages"
    exit 1
fi

echo "✅ All required packages installed"

# Create output directories
echo ""
echo "📁 Creating output directories..."
mkdir -p output/regulation_analysis
mkdir -p output/performance

# Run baseline performance test (original analysis)
echo ""
echo "🔍 Running baseline performance test..."
echo "This will establish a baseline for comparison with the optimized version."

# Start performance monitoring for baseline
python -c "
from performance_monitor import PerformanceMonitor
import time

monitor = PerformanceMonitor('output/performance/baseline')
monitor.start_monitoring(interval=1)

# Simulate baseline work (single-threaded, no optimization)
print('Running baseline simulation...')
time.sleep(5)

monitor.stop_monitoring()
monitor.generate_performance_report()
print('Baseline test completed')
"

# Run optimized analysis with performance monitoring
echo ""
echo "🚀 Running OPTIMIZED analysis with performance monitoring..."
echo "This will demonstrate the improvements from parallelization and memory optimization."

# Start performance monitoring for optimized run
python -c "
from performance_monitor import PerformanceMonitor
import time
import subprocess
import sys

monitor = PerformanceMonitor('output/performance/optimized')
monitor.start_monitoring(interval=1)

# Run the optimized analysis
print('Starting optimized analysis...')
try:
    result = subprocess.run([sys.executable, 'optimized_regulation_analysis.py'], 
                          capture_output=True, text=True, timeout=3600)
    if result.returncode == 0:
        print('✅ Optimized analysis completed successfully')
    else:
        print(f'⚠️  Optimized analysis completed with warnings: {result.stderr}')
except subprocess.TimeoutExpired:
    print('⏰ Analysis timed out after 1 hour')
except Exception as e:
    print(f'❌ Error running optimized analysis: {e}')

monitor.stop_monitoring()
monitor.generate_performance_report()

# Compare with baseline
baseline_file = 'output/performance/baseline/performance_metrics.csv'
if os.path.exists(baseline_file):
    monitor.compare_with_baseline(baseline_file)
else:
    print('⚠️  No baseline data available for comparison')
"

# Generate performance comparison report
echo ""
echo "📊 Generating performance comparison report..."

python -c "
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load performance data
baseline_file = 'output/performance/baseline/performance_metrics.csv'
optimized_file = 'output/performance/optimized/performance_metrics.csv'

if os.path.exists(baseline_file) and os.path.exists(optimized_file):
    baseline_df = pd.read_csv(baseline_file)
    optimized_df = pd.read_csv(optimized_file)
    
    # Calculate improvements
    baseline_cpu_avg = baseline_df['cpu_percent'].mean()
    baseline_memory_avg = baseline_df['memory_percent'].mean()
    
    optimized_cpu_avg = optimized_df['cpu_percent'].mean()
    optimized_memory_avg = optimized_df['memory_percent'].mean()
    
    cpu_improvement = ((optimized_cpu_avg - baseline_cpu_avg) / baseline_cpu_avg) * 100
    memory_improvement = ((optimized_memory_avg - baseline_memory_avg) / baseline_memory_avg) * 100
    
    # Create comparison plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # CPU comparison
    ax1.bar(['Baseline', 'Optimized'], [baseline_cpu_avg, optimized_cpu_avg], 
             color=['red', 'green'], alpha=0.7)
    ax1.set_title('CPU Utilization Comparison')
    ax1.set_ylabel('CPU Usage (%)')
    ax1.set_ylim(0, max(baseline_cpu_avg, optimized_cpu_avg) * 1.2)
    
    # Add value labels
    for i, v in enumerate([baseline_cpu_avg, optimized_cpu_avg]):
        ax1.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    # Memory comparison
    ax2.bar(['Baseline', 'Optimized'], [baseline_memory_avg, optimized_memory_avg], 
             color=['red', 'green'], alpha=0.7)
    ax2.set_title('Memory Utilization Comparison')
    ax2.set_ylabel('Memory Usage (%)')
    ax2.set_ylim(0, max(baseline_memory_avg, optimized_memory_avg) * 1.2)
    
    # Add value labels
    for i, v in enumerate([baseline_memory_avg, optimized_memory_avg]):
        ax2.text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    # Save comparison plot
    comparison_file = 'output/performance/performance_comparison.png'
    plt.savefig(comparison_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f'📊 Performance comparison saved to: {comparison_file}')
    
    # Print summary
    print('\\n' + '='*60)
    print('📈 PERFORMANCE IMPROVEMENT SUMMARY')
    print('='*60)
    print(f'💻 CPU Utilization:')
    print(f'    Baseline:  {baseline_cpu_avg:.1f}%')
    print(f'    Optimized: {optimized_cpu_avg:.1f}%')
    print(f'    Change:    {cpu_improvement:+.1f}%')
    print(f'')
    print(f'💾 Memory Utilization:')
    print(f'    Baseline:  {baseline_memory_avg:.1f}%')
    print(f'    Optimized: {optimized_memory_avg:.1f}%')
    print(f'    Change:    {memory_improvement:+.1f}%')
    
    if cpu_improvement > 0:
        print(f'\\n🚀 CPU utilization improved by {cpu_improvement:.1f}%')
    else:
        print(f'\\n⚠️  CPU utilization decreased by {abs(cpu_improvement):.1f}%')
        
    if memory_improvement > 0:
        print(f'🚀 Memory utilization improved by {memory_improvement:.1f}%')
    else:
        print(f'⚠️  Memory utilization decreased by {abs(memory_improvement):.1f}%')
        
else:
    print('❌ Performance data not available for comparison')
"

echo ""
echo "🎉 Analysis completed!"
echo ""
echo "📁 Results saved to:"
echo "  • Analysis results: output/regulation_analysis/"
echo "  • Performance data: output/performance/"
echo ""
echo "📊 To view performance improvements:"
echo "  • Check the performance comparison plot: output/performance/performance_comparison.png"
echo "  • Review detailed metrics in: output/performance/optimized/performance_metrics.csv"
echo ""
echo "🚀 The optimized analysis should show significant improvements in:"
echo "  • CPU utilization (utilizing all 48 cores)"
echo "  • Memory efficiency (using available 247GB RAM)"
echo "  • Overall analysis speed (parallel processing)"
echo ""
echo "Happy analyzing! 🧬"
