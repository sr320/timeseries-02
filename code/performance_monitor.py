#!/usr/bin/env python3
"""
Performance Monitor for Timeseries Analysis

This script monitors system resources during analysis execution to demonstrate
the improvements from parallelization and memory optimization.
"""

import psutil
import time
import threading
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import os

class PerformanceMonitor:
    def __init__(self, output_dir="output/performance"):
        """Initialize the performance monitor."""
        self.output_dir = output_dir
        self.monitoring = False
        self.metrics = {
            'timestamps': [],
            'cpu_percent': [],
            'memory_percent': [],
            'memory_used_gb': [],
            'disk_io_read': [],
            'disk_io_write': [],
            'network_sent': [],
            'network_recv': []
        }
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Get initial system info
        self.system_info = self._get_system_info()
        
    def _get_system_info(self):
        """Get system information."""
        info = {
            'cpu_count': psutil.cpu_count(),
            'cpu_freq': psutil.cpu_freq(),
            'memory_total_gb': psutil.virtual_memory().total / (1024**3),
            'disk_usage': psutil.disk_usage('/'),
            'python_version': f"{psutil.sys.version_info.major}.{psutil.sys.version_info.minor}.{psutil.sys.version_info.micro}"
        }
        return info
        
    def start_monitoring(self, interval=1):
        """Start monitoring system resources."""
        self.monitoring = True
        self.monitor_thread = threading.Thread(target=self._monitor_loop, args=(interval,))
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
        print(f"🚀 Performance monitoring started (interval: {interval}s)")
        print(f"💻 System: {self.system_info['cpu_count']} CPU cores, {self.system_info['memory_total_gb']:.1f} GB RAM")
        
    def stop_monitoring(self):
        """Stop monitoring system resources."""
        self.monitoring = False
        if hasattr(self, 'monitor_thread'):
            self.monitor_thread.join()
        print("⏹️  Performance monitoring stopped")
        
    def _monitor_loop(self, interval):
        """Main monitoring loop."""
        while self.monitoring:
            try:
                # CPU usage
                cpu_percent = psutil.cpu_percent(interval=interval, percpu=False)
                
                # Memory usage
                memory = psutil.virtual_memory()
                memory_percent = memory.percent
                memory_used_gb = memory.used / (1024**3)
                
                # Disk I/O
                disk_io = psutil.disk_io_counters()
                disk_read = disk_io.read_bytes if disk_io else 0
                disk_write = disk_io.write_bytes if disk_io else 0
                
                # Network I/O
                network_io = psutil.net_io_counters()
                network_sent = network_io.bytes_sent if network_io else 0
                network_recv = network_io.bytes_recv if network_io else 0
                
                # Store metrics
                timestamp = time.time()
                self.metrics['timestamps'].append(timestamp)
                self.metrics['cpu_percent'].append(cpu_percent)
                self.metrics['memory_percent'].append(memory_percent)
                self.metrics['memory_used_gb'].append(memory_used_gb)
                self.metrics['disk_io_read'].append(disk_read)
                self.metrics['disk_io_write'].append(disk_write)
                self.metrics['network_sent'].append(network_sent)
                self.metrics['network_recv'].append(network_recv)
                
                # Print current status
                print(f"📊 CPU: {cpu_percent:5.1f}% | RAM: {memory_percent:5.1f}% ({memory_used_gb:6.1f} GB) | Time: {datetime.fromtimestamp(timestamp).strftime('%H:%M:%S')}")
                
            except Exception as e:
                print(f"❌ Monitoring error: {e}")
                time.sleep(interval)
                
    def get_current_stats(self):
        """Get current system statistics."""
        if not self.metrics['timestamps']:
            return None
            
        return {
            'cpu_percent': self.metrics['cpu_percent'][-1] if self.metrics['cpu_percent'] else 0,
            'memory_percent': self.metrics['memory_percent'][-1] if self.metrics['memory_percent'] else 0,
            'memory_used_gb': self.metrics['memory_used_gb'][-1] if self.metrics['memory_used_gb'] else 0,
            'monitoring_duration': self.metrics['timestamps'][-1] - self.metrics['timestamps'][0] if len(self.metrics['timestamps']) > 1 else 0
        }
        
    def generate_performance_report(self):
        """Generate comprehensive performance report."""
        if not self.metrics['timestamps']:
            print("❌ No performance data available")
            return
            
        print("\n" + "="*60)
        print("📊 PERFORMANCE MONITORING REPORT")
        print("="*60)
        
        # Calculate statistics
        cpu_avg = np.mean(self.metrics['cpu_percent'])
        cpu_max = np.max(self.metrics['cpu_percent'])
        memory_avg = np.mean(self.metrics['memory_percent'])
        memory_max = np.max(self.metrics['memory_percent'])
        memory_used_max = np.max(self.metrics['memory_used_gb'])
        
        duration = self.metrics['timestamps'][-1] - self.metrics['timestamps'][0]
        
        print(f"⏱️  Monitoring Duration: {duration:.1f} seconds ({duration/60:.1f} minutes)")
        print(f"💻 CPU Usage:")
        print(f"    Average: {cpu_avg:.1f}%")
        print(f"    Maximum: {cpu_max:.1f}%")
        print(f"    Cores Utilized: {cpu_max/100 * self.system_info['cpu_count']:.1f}/{self.system_info['cpu_count']}")
        print(f"💾 Memory Usage:")
        print(f"    Average: {memory_avg:.1f}%")
        print(f"    Maximum: {memory_max:.1f}%")
        print(f"    Peak Usage: {memory_used_max:.1f} GB / {self.system_info['memory_total_gb']:.1f} GB")
        print(f"    Memory Efficiency: {memory_used_max/self.system_info['memory_total_gb']*100:.1f}%")
        
        # Resource utilization assessment
        print(f"\n🎯 RESOURCE UTILIZATION ASSESSMENT:")
        if cpu_max > 80:
            print(f"    ✅ CPU: Excellent utilization ({cpu_max:.1f}% peak)")
        elif cpu_max > 50:
            print(f"    ⚠️  CPU: Good utilization ({cpu_max:.1f}% peak)")
        else:
            print(f"    ❌ CPU: Low utilization ({cpu_max:.1f}% peak) - consider increasing parallelization")
            
        if memory_max > 80:
            print(f"    ✅ Memory: Excellent utilization ({memory_max:.1f}% peak)")
        elif memory_max > 50:
            print(f"    ⚠️  Memory: Good utilization ({memory_max:.1f}% peak)")
        else:
            print(f"    ❌ Memory: Low utilization ({memory_max:.1f}% peak) - consider larger batch sizes")
            
        # Save detailed metrics
        self._save_metrics()
        
    def _save_metrics(self):
        """Save detailed metrics to files."""
        # Save raw metrics
        metrics_file = os.path.join(self.output_dir, "performance_metrics.csv")
        import pandas as pd
        
        df = pd.DataFrame({
            'timestamp': [datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S') for ts in self.metrics['timestamps']],
            'cpu_percent': self.metrics['cpu_percent'],
            'memory_percent': self.metrics['memory_percent'],
            'memory_used_gb': self.metrics['memory_used_gb'],
            'disk_io_read_mb': [x/(1024**2) for x in self.metrics['disk_io_read']],
            'disk_io_write_mb': [x/(1024**2) for x in self.metrics['disk_io_write']],
            'network_sent_mb': [x/(1024**2) for x in self.metrics['network_sent']],
            'network_recv_mb': [x/(1024**2) for x in self.metrics['network_recv']]
        })
        
        df.to_csv(metrics_file, index=False)
        print(f"📁 Detailed metrics saved to: {metrics_file}")
        
        # Generate performance plots
        self._generate_performance_plots()
        
    def _generate_performance_plots(self):
        """Generate performance visualization plots."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle('System Performance During Analysis', fontsize=16)
        
        # CPU Usage
        axes[0, 0].plot(self.metrics['timestamps'], self.metrics['cpu_percent'], 'b-', linewidth=2)
        axes[0, 0].set_title('CPU Usage Over Time')
        axes[0, 0].set_ylabel('CPU Usage (%)')
        axes[0, 0].set_ylim(0, 100)
        axes[0, 0].grid(True, alpha=0.3)
        
        # Memory Usage
        axes[0, 1].plot(self.metrics['timestamps'], self.metrics['memory_percent'], 'r-', linewidth=2)
        axes[0, 1].set_title('Memory Usage Over Time')
        axes[0, 1].set_ylabel('Memory Usage (%)')
        axes[0, 1].set_ylim(0, 100)
        axes[0, 1].grid(True, alpha=0.3)
        
        # Memory Used in GB
        axes[1, 0].plot(self.metrics['timestamps'], self.metrics['memory_used_gb'], 'g-', linewidth=2)
        axes[1, 0].set_title('Memory Usage (GB)')
        axes[1, 0].set_ylabel('Memory Used (GB)')
        axes[1, 0].axhline(y=self.system_info['memory_total_gb'], color='r', linestyle='--', alpha=0.7, label='Total RAM')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Resource Utilization Summary
        cpu_avg = np.mean(self.metrics['cpu_percent'])
        memory_avg = np.mean(self.metrics['memory_percent'])
        
        axes[1, 1].bar(['CPU', 'Memory'], [cpu_avg, memory_avg], color=['blue', 'red'], alpha=0.7)
        axes[1, 1].set_title('Average Resource Utilization')
        axes[1, 1].set_ylabel('Utilization (%)')
        axes[1, 1].set_ylim(0, 100)
        axes[1, 1].grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, v in enumerate([cpu_avg, memory_avg]):
            axes[1, 1].text(i, v + 1, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        plot_file = os.path.join(self.output_dir, "performance_analysis.png")
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 Performance plots saved to: {plot_file}")
        
    def compare_with_baseline(self, baseline_file):
        """Compare current performance with a baseline run."""
        if not os.path.exists(baseline_file):
            print(f"❌ Baseline file not found: {baseline_file}")
            return
            
        try:
            import pandas as pd
            baseline_df = pd.read_csv(baseline_file)
            
            print("\n" + "="*60)
            print("📈 PERFORMANCE COMPARISON WITH BASELINE")
            print("="*60)
            
            # Calculate current averages
            current_cpu_avg = np.mean(self.metrics['cpu_percent'])
            current_memory_avg = np.mean(self.metrics['memory_percent'])
            
            # Calculate baseline averages
            baseline_cpu_avg = baseline_df['cpu_percent'].mean()
            baseline_memory_avg = baseline_df['memory_percent'].mean()
            
            # Calculate improvements
            cpu_improvement = ((current_cpu_avg - baseline_cpu_avg) / baseline_cpu_avg) * 100
            memory_improvement = ((current_memory_avg - baseline_memory_avg) / baseline_memory_avg) * 100
            
            print(f"💻 CPU Utilization:")
            print(f"    Baseline: {baseline_cpu_avg:.1f}%")
            print(f"    Current:  {current_cpu_avg:.1f}%")
            print(f"    Change:   {cpu_improvement:+.1f}%")
            
            print(f"💾 Memory Utilization:")
            print(f"    Baseline: {baseline_memory_avg:.1f}%")
            print(f"    Current:  {current_memory_avg:.1f}%")
            print(f"    Change:   {memory_improvement:+.1f}%")
            
            if cpu_improvement > 0:
                print(f"🚀 CPU utilization improved by {cpu_improvement:.1f}%")
            else:
                print(f"⚠️  CPU utilization decreased by {abs(cpu_improvement):.1f}%")
                
            if memory_improvement > 0:
                print(f"🚀 Memory utilization improved by {memory_improvement:.1f}%")
            else:
                print(f"⚠️  Memory utilization decreased by {abs(memory_improvement):.1f}%")
                
        except Exception as e:
            print(f"❌ Error comparing with baseline: {e}")

def main():
    """Example usage of the performance monitor."""
    print("🔍 Performance Monitor for Timeseries Analysis")
    print("="*50)
    
    # Initialize monitor
    monitor = PerformanceMonitor()
    
    # Start monitoring
    monitor.start_monitoring(interval=2)
    
    try:
        # Simulate some work
        print("\n🔄 Simulating analysis work...")
        time.sleep(10)  # Monitor for 10 seconds
        
    except KeyboardInterrupt:
        print("\n⏹️  Monitoring interrupted by user")
    finally:
        # Stop monitoring
        monitor.stop_monitoring()
        
        # Generate report
        monitor.generate_performance_report()

if __name__ == "__main__":
    main()
