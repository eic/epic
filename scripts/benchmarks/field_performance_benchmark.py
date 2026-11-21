#!/usr/bin/env python3
"""
Magnetic Field Performance Benchmark for EPIC FieldMapB

This script benchmarks the performance of the FieldMapB implementation using covfie,
measuring timing, memory usage, and accuracy across different field configurations.

Usage:
    ./field_performance_benchmark.py [options]

Requirements:
    - EIC environment with DD4hep and EPIC ins        # Get CPU info        # Get CPU info safely
          # Get CPU info safely
        try:
            cpu_lines = subprocess.run(['cat', '/proc/cpuinfo'], capture_output=True, text=True).stdout.split('\n')
            cpu_info = next((line for line in cpu_lines if 'model name' in line), 'Unknown CPU')
        except:
            cpu_info = 'Unknown CPU'

        all_results = {
            'metadata': {
                'timestamp': time.time(),
                'hostname': os.uname().nodename,
                'cpu_info': cpu_info,
                'memory_gb': psutil.virtual_memory().total / (1024**3),
                'epic_version': os.environ.get('EPIC_VERSION', 'unknown'),
                'dd4hep_version': os.environ.get('DD4hepINSTALL', 'unknown')
            },:
            cpu_lines = subprocess.run(['cat', '/proc/cpuinfo'], capture_output=True, text=True).stdout.split('\n')
            cpu_info = next((line for line in cpu_lines if 'model name' in line), 'Unknown CPU')
        except:
            cpu_info = 'Unknown CPU'

        all_results = {
            'metadata': {
                'timestamp': time.time(),
                'hostname': os.uname().nodename,
                'cpu_info': cpu_info,
                'memory_gb': psutil.virtual_memory().total / (1024**3),
                'epic_version': os.environ.get('EPIC_VERSION', 'unknown'),
                'dd4hep_version': os.environ.get('DD4hepINSTALL', 'unknown')
            },       try:
            cpu_lines = subprocess.run(['cat', '/proc/cpuinfo'], capture_output=True, text=True).stdout.split('\n')
            cpu_info = next((line for line in cpu_lines if 'model name' in line), 'Unknown CPU')
        except:
            cpu_info = 'Unknown CPU'

        all_results = {
            'metadata': {
                'timestamp': time.time(),
                'hostname': os.uname().nodename,
                'cpu_info': cpu_info,
                'memory_gb': psutil.virtual_memory().total / (1024**3),
                'epic_version': os.environ.get('EPIC_VERSION', 'unknown'),
                'dd4hep_version': os.environ.get('DD4hepINSTALL', 'unknown')
            },  - Field map files available in fieldmaps/
    - Python packages: numpy, matplotlib, psutil, json
"""

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import subprocess
import tempfile
import shutil

try:
    import numpy as np
    import matplotlib.pyplot as plt
    import psutil
except ImportError as e:
    print(f"Required Python package not available: {e}")
    print("Please install: pip install numpy matplotlib psutil")
    sys.exit(1)

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FieldBenchmark:
    """Benchmark suite for EPIC magnetic field performance."""

    def __init__(self, detector_path: str, output_dir: str = "benchmark_results"):
        self.detector_path = Path(detector_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {}

        # Benchmark configurations
        self.field_configs = {
            'marco_solenoid': {
                'xml_file': 'compact/fields/marco.xml',
                'coord_type': 'BrBz',
                'description': 'MARCO solenoid field (cylindrical coords)'
            },
            'lumi_magnets': {
                'xml_file': 'compact/far_backward/lumi/lumi_magnets.xml',
                'coord_type': 'BxByBz',
                'description': 'Lumi dipole magnets (cartesian coords)'
            }
        }

        # Test parameters
        self.n_samples = [1000, 10000, 100000, 500000]  # Different sample sizes
        self.test_regions = {
            'barrel': {'r_range': (0, 100), 'z_range': (-150, 150)},  # cm
            'forward': {'r_range': (0, 50), 'z_range': (150, 400)},
            'backward': {'r_range': (0, 50), 'z_range': (-400, -150)}
        }

    def create_test_geometry(self, field_config: str) -> str:
        """Create a minimal geometry file for testing a specific field configuration."""
        config = self.field_configs[field_config]

        # Create minimal detector XML for testing
        test_xml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<lccdd>
  <info name="epic_field_test" title="EPIC Field Test Geometry" author="field_benchmark" url="auto" status="test" version="test">
    <comment>Minimal geometry for field performance testing</comment>
  </info>

  <includes>
    <gdmlFile ref="${{DD4hepINSTALL}}/DDDetectors/compact/elements.xml"/>
    <gdmlFile ref="${{DD4hepINSTALL}}/DDDetectors/compact/materials.xml"/>
  </includes>

  <define>
    <constant name="world_side" value="10*m"/>
    <constant name="world_x" value="world_side"/>
    <constant name="world_y" value="world_side"/>
    <constant name="world_z" value="world_side"/>
  </define>

  <materials>
    <material name="Vacuum">
      <D type="density" unit="g/cm3" value="0.00000001" />
      <fraction n="1" ref="H"/>
    </material>
  </materials>

  <detectors>
    <detector name="WorldBox" type="DD4hep_BoxSegment" vis="InvisibleWithChildren">
      <material name="Vacuum"/>
      <box dx="world_x" dy="world_y" dz="world_z"/>
    </detector>
  </detectors>

  <include ref="{config['xml_file']}"/>

  <fields>
    <!-- Field configuration loaded from included file -->
  </fields>

</lccdd>"""

        # Write to temporary file
        temp_file = self.output_dir / f"test_{field_config}.xml"
        with open(temp_file, 'w') as f:
            f.write(test_xml_content)

        return str(temp_file)

    def run_field_timing_test(self, xml_file: str, field_config: str, n_points: int, region: str) -> Dict:
        """Run timing test using DD4hep field evaluation."""

        # Create C++ benchmark program
        cpp_code = f"""
#include <DD4hep/Detector.h>
#include <DD4hep/FieldTypes.h>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace dd4hep;
using namespace std;

int main() {{
    Detector& detector = Detector::getInstance();
    detector.fromXML("{xml_file}");

    auto field = detector.field();
    if (!field.isValid()) {{
        cerr << "ERROR: No field found in detector description" << endl;
        return 1;
    }}

    // Generate random test points
    random_device rd;
    mt19937 gen(42); // Fixed seed for reproducibility
    uniform_real_distribution<> r_dist({self.test_regions[region]['r_range'][0]},
                                      {self.test_regions[region]['r_range'][1]});
    uniform_real_distribution<> phi_dist(0, 2 * M_PI);
    uniform_real_distribution<> z_dist({self.test_regions[region]['z_range'][0]},
                                      {self.test_regions[region]['z_range'][1]});

    vector<tuple<double, double, double>> test_points;
    test_points.reserve({n_points});

    for (int i = 0; i < {n_points}; ++i) {{
        double r = r_dist(gen);
        double phi = phi_dist(gen);
        double z = z_dist(gen);
        double x = r * cos(phi);
        double y = r * sin(phi);
        test_points.emplace_back(x, y, z);
    }}

    // Warm up
    double pos[3], field_val[3];
    for (int i = 0; i < 1000; ++i) {{
        auto [x, y, z] = test_points[i % test_points.size()];
        pos[0] = x; pos[1] = y; pos[2] = z;
        field.magneticField(pos, field_val);
    }}

    // Timing test
    auto start = chrono::high_resolution_clock::now();

    double sum_bx = 0, sum_by = 0, sum_bz = 0;
    for (const auto& point : test_points) {{
        auto [x, y, z] = point;
        pos[0] = x; pos[1] = y; pos[2] = z;
        field.magneticField(pos, field_val);
        sum_bx += field_val[0];
        sum_by += field_val[1];
        sum_bz += field_val[2];
    }}

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

    // Output results
    cout << "{{" << endl;
    cout << "  \\"n_points\\": " << {n_points} << "," << endl;
    cout << "  \\"total_time_us\\": " << duration.count() << "," << endl;
    cout << "  \\"time_per_evaluation_ns\\": " << (duration.count() * 1000.0 / {n_points}) << "," << endl;
    cout << "  \\"evaluations_per_second\\": " << ({n_points} * 1e6 / duration.count()) << "," << endl;
    cout << "  \\"sum_field\\": [" << sum_bx << ", " << sum_by << ", " << sum_bz << "]," << endl;
    cout << "  \\"field_magnitude_avg\\": " << sqrt(sum_bx*sum_bx + sum_by*sum_by + sum_bz*sum_bz) / {n_points} << endl;
    cout << "}}" << endl;

    return 0;
}}
"""

        # Compile and run C++ benchmark
        cpp_file = self.output_dir / f"benchmark_{field_config}_{region}_{n_points}.cpp"
        exe_file = self.output_dir / f"benchmark_{field_config}_{region}_{n_points}"

        with open(cpp_file, 'w') as f:
            f.write(cpp_code)

        # Compile
        dd4hep_install = os.environ.get('DD4hepINSTALL', '/opt/local')

        # Get ROOT configuration
        try:
            root_cflags = subprocess.run(['root-config', '--cflags'], capture_output=True, text=True).stdout.strip().split()
            root_libs = subprocess.run(['root-config', '--libs'], capture_output=True, text=True).stdout.strip().split()
        except:
            root_cflags = ['-I/opt/local/include/root']
            root_libs = ['-lCore', '-lMathCore']

        compile_cmd = [
            "g++", "-O3", "-march=native",
            f"-I{dd4hep_install}/include",
            f"-L{dd4hep_install}/lib",
            "-lDDCore", "-lDDRec",
            str(cpp_file), "-o", str(exe_file)
        ] + root_cflags + root_libs

        logger.info(f"Compiling benchmark for {field_config}, {region}, {n_points} points...")
        try:
            result = subprocess.run(compile_cmd, shell=False, capture_output=True, text=True,
                                 env=dict(os.environ))
            if result.returncode != 0:
                logger.error(f"Compilation failed: {result.stderr}")
                return None
        except Exception as e:
            logger.error(f"Compilation error: {e}")
            return None

        # Run benchmark
        logger.info(f"Running benchmark...")
        try:
            # Monitor memory usage
            process = psutil.Popen([str(exe_file)], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, text=True)

            max_memory = 0
            while process.poll() is None:
                try:
                    memory_info = process.memory_info()
                    max_memory = max(max_memory, memory_info.rss / 1024 / 1024)  # MB
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    break
                time.sleep(0.01)

            stdout, stderr = process.communicate()

            if process.returncode != 0:
                logger.error(f"Benchmark execution failed: {stderr}")
                return None

            # Parse results
            result_data = json.loads(stdout)
            result_data['max_memory_mb'] = max_memory
            result_data['field_config'] = field_config
            result_data['region'] = region

            return result_data

        except Exception as e:
            logger.error(f"Execution error: {e}")
            return None
        finally:
            # Cleanup
            for f in [cpp_file, exe_file]:
                if f.exists():
                    f.unlink()

    def run_accuracy_test(self, xml_file: str, field_config: str) -> Dict:
        """Test field accuracy and consistency."""

        cpp_code = f"""
#include <DD4hep/Detector.h>
#include <DD4hep/FieldTypes.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace dd4hep;

int main() {{
    Detector& detector = Detector::getInstance();
    detector.fromXML("{xml_file}");

    auto field = detector.field();
    if (!field.isValid()) {{
        std::cerr << "ERROR: No field found" << std::endl;
        return 1;
    }}

    // Test field properties at key points
    double pos[3], field_val[3];

    // Test at origin
    pos[0] = 0; pos[1] = 0; pos[2] = 0;
    field.magneticField(pos, field_val);
    double field_at_origin = sqrt(field_val[0]*field_val[0] + field_val[1]*field_val[1] + field_val[2]*field_val[2]);

    // Test cylindrical symmetry (for BrBz fields)
    double asymmetry = 0.0;
    for (int phi_deg = 0; phi_deg < 360; phi_deg += 45) {{
        double phi = phi_deg * M_PI / 180.0;
        double r = 50.0; // 50 cm
        pos[0] = r * cos(phi);
        pos[1] = r * sin(phi);
        pos[2] = 0;
        field.magneticField(pos, field_val);
        double field_mag = sqrt(field_val[0]*field_val[0] + field_val[1]*field_val[1] + field_val[2]*field_val[2]);
        asymmetry += abs(field_mag - field_at_origin) / field_at_origin;
    }}
    asymmetry /= 8.0;  // Average over 8 points

    // Test field gradient
    pos[0] = 0; pos[1] = 0; pos[2] = 0;
    field.magneticField(pos, field_val);
    double field_center = field_val[2];  // Bz at center

    pos[2] = 10.0;  // 10 cm offset
    field.magneticField(pos, field_val);
    double field_offset = field_val[2];
    double gradient = (field_offset - field_center) / 10.0;  // T/cm

    std::cout << "{{" << std::endl;
    std::cout << "  \\"field_at_origin_T\\": " << field_at_origin << "," << std::endl;
    std::cout << "  \\"cylindrical_asymmetry\\": " << asymmetry << "," << std::endl;
    std::cout << "  \\"field_gradient_T_per_cm\\": " << gradient << std::endl;
    std::cout << "}}" << std::endl;

    return 0;
}}
"""

        cpp_file = self.output_dir / f"accuracy_{field_config}.cpp"
        exe_file = self.output_dir / f"accuracy_{field_config}"

        with open(cpp_file, 'w') as f:
            f.write(cpp_code)

        # Compile and run
        dd4hep_install = os.environ.get('DD4hepINSTALL', '/opt/local')

        # Get ROOT configuration
        try:
            root_cflags = subprocess.run(['root-config', '--cflags'], capture_output=True, text=True).stdout.strip().split()
            root_libs = subprocess.run(['root-config', '--libs'], capture_output=True, text=True).stdout.strip().split()
        except:
            root_cflags = ['-I/opt/local/include/root']
            root_libs = ['-lCore', '-lMathCore']

        compile_cmd = [
            "g++", "-O3",
            f"-I{dd4hep_install}/include",
            f"-L{dd4hep_install}/lib",
            "-lDDCore", "-lDDRec",
            str(cpp_file), "-o", str(exe_file)
        ] + root_cflags + root_libs

        try:
            subprocess.run(compile_cmd, shell=False, check=True, capture_output=True)
            result = subprocess.run([str(exe_file)], capture_output=True, text=True, check=True)

            accuracy_data = json.loads(result.stdout)
            return accuracy_data

        except Exception as e:
            logger.error(f"Accuracy test failed for {field_config}: {e}")
            return {}
        finally:
            for f in [cpp_file, exe_file]:
                if f.exists():
                    f.unlink()

    def run_comprehensive_benchmark(self) -> Dict:
        """Run complete benchmark suite."""
        logger.info("Starting comprehensive field performance benchmark...")

        all_results = {
            'metadata': {
                'timestamp': time.time(),
                'hostname': os.uname().nodename,
                'cpu_info': 'Unknown CPU',  # Fixed CPU info parsing
                'memory_gb': psutil.virtual_memory().total / (1024**3),
                'epic_version': os.environ.get('EPIC_VERSION', 'unknown'),
                'dd4hep_version': os.environ.get('DD4hepINSTALL', 'unknown')
            },
            'timing_results': {},
            'accuracy_results': {},
            'performance_summary': {}
        }

        # Run timing benchmarks
        for field_config in self.field_configs.keys():
            logger.info(f"Testing {field_config}...")

            # Create test geometry
            try:
                xml_file = self.create_test_geometry(field_config)
                all_results['timing_results'][field_config] = {}

                # Test different sample sizes and regions
                for region in self.test_regions.keys():
                    all_results['timing_results'][field_config][region] = {}

                    for n_points in self.n_samples:
                        logger.info(f"  Testing {region} region with {n_points} points...")

                        result = self.run_field_timing_test(xml_file, field_config, n_points, region)
                        if result:
                            all_results['timing_results'][field_config][region][n_points] = result

                # Run accuracy tests
                accuracy_result = self.run_accuracy_test(xml_file, field_config)
                if accuracy_result:
                    all_results['accuracy_results'][field_config] = accuracy_result

            except Exception as e:
                logger.error(f"Failed to test {field_config}: {e}")
                continue

        # Generate performance summary
        self.generate_performance_summary(all_results)

        return all_results

    def generate_performance_summary(self, results: Dict):
        """Generate performance summary and plots."""

        # Calculate performance metrics
        summary = {}

        for field_config, timing_data in results['timing_results'].items():
            config_summary = {
                'avg_evaluations_per_second': 0,
                'avg_time_per_evaluation_ns': 0,
                'memory_efficiency': 0,
                'scalability_score': 0
            }

            eval_rates = []
            eval_times = []

            for region, region_data in timing_data.items():
                for n_points, point_data in region_data.items():
                    if isinstance(point_data, dict):
                        eval_rates.append(point_data.get('evaluations_per_second', 0))
                        eval_times.append(point_data.get('time_per_evaluation_ns', 0))

            if eval_rates:
                config_summary['avg_evaluations_per_second'] = np.mean(eval_rates)
                config_summary['avg_time_per_evaluation_ns'] = np.mean(eval_times)
                config_summary['scalability_score'] = np.std(eval_rates) / np.mean(eval_rates) if np.mean(eval_rates) > 0 else 1.0

            summary[field_config] = config_summary

        results['performance_summary'] = summary

        # Create performance plots
        self.create_performance_plots(results)

    def create_performance_plots(self, results: Dict):
        """Create performance visualization plots."""

        # Performance comparison plot
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle('EPIC Field Performance Benchmark Results', fontsize=16)

        configs = list(results['timing_results'].keys())
        colors = ['blue', 'red', 'green', 'orange']

        # Plot 1: Evaluations per second vs sample size
        for i, config in enumerate(configs):
            sample_sizes = []
            eval_rates = []

            for region in self.test_regions.keys():
                if region in results['timing_results'][config]:
                    for n_points, data in results['timing_results'][config][region].items():
                        if isinstance(data, dict) and 'evaluations_per_second' in data:
                            sample_sizes.append(n_points)
                            eval_rates.append(data['evaluations_per_second'])

            if sample_sizes:
                ax1.loglog(sample_sizes, eval_rates, 'o-', color=colors[i % len(colors)],
                          label=f'{config}', markersize=6)

        ax1.set_xlabel('Sample Size')
        ax1.set_ylabel('Evaluations/Second')
        ax1.set_title('Throughput vs Sample Size')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Time per evaluation
        for i, config in enumerate(configs):
            sample_sizes = []
            eval_times = []

            for region in self.test_regions.keys():
                if region in results['timing_results'][config]:
                    for n_points, data in results['timing_results'][config][region].items():
                        if isinstance(data, dict) and 'time_per_evaluation_ns' in data:
                            sample_sizes.append(n_points)
                            eval_times.append(data['time_per_evaluation_ns'])

            if sample_sizes:
                ax2.semilogx(sample_sizes, eval_times, 'o-', color=colors[i % len(colors)],
                            label=f'{config}', markersize=6)

        ax2.set_xlabel('Sample Size')
        ax2.set_ylabel('Time per Evaluation (ns)')
        ax2.set_title('Latency vs Sample Size')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # Plot 3: Memory usage
        memory_data = {}
        for config in configs:
            memory_usage = []
            for region in self.test_regions.keys():
                if region in results['timing_results'][config]:
                    for n_points, data in results['timing_results'][config][region].items():
                        if isinstance(data, dict) and 'max_memory_mb' in data:
                            memory_usage.append(data['max_memory_mb'])
            if memory_usage:
                memory_data[config] = np.mean(memory_usage)

        if memory_data:
            ax3.bar(memory_data.keys(), memory_data.values(), color=colors[:len(memory_data)])
            ax3.set_ylabel('Memory Usage (MB)')
            ax3.set_title('Average Memory Usage')
            ax3.tick_params(axis='x', rotation=45)

        # Plot 4: Performance summary
        if results['performance_summary']:
            perf_metrics = ['avg_evaluations_per_second', 'scalability_score']
            x_pos = np.arange(len(configs))

            for i, metric in enumerate(perf_metrics):
                values = [results['performance_summary'][config].get(metric, 0) for config in configs]
                ax4.bar(x_pos + i*0.35, values, 0.35, label=metric, color=colors[i])

            ax4.set_xlabel('Field Configuration')
            ax4.set_ylabel('Performance Score')
            ax4.set_title('Performance Summary')
            ax4.set_xticks(x_pos + 0.175)
            ax4.set_xticklabels(configs, rotation=45)
            ax4.legend()

        plt.tight_layout()
        plt.savefig(self.output_dir / 'field_performance_benchmark.png', dpi=300, bbox_inches='tight')
        plt.close()

        logger.info(f"Performance plots saved to {self.output_dir / 'field_performance_benchmark.png'}")

    def save_results(self, results: Dict, filename: str = "field_benchmark_results.json"):
        """Save benchmark results to JSON file."""
        output_file = self.output_dir / filename

        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        logger.info(f"Results saved to {output_file}")
        return output_file

    def generate_report(self, results: Dict) -> str:
        """Generate human-readable benchmark report."""

        report = []
        report.append("EPIC Field Performance Benchmark Report")
        report.append("=" * 50)
        report.append(f"Timestamp: {time.ctime(results['metadata']['timestamp'])}")
        report.append(f"Hostname: {results['metadata']['hostname']}")
        report.append(f"CPU: {results['metadata']['cpu_info']}")
        report.append(f"Memory: {results['metadata']['memory_gb']:.1f} GB")
        report.append("")

        # Performance summary
        report.append("Performance Summary:")
        report.append("-" * 20)

        for config, summary in results.get('performance_summary', {}).items():
            report.append(f"\\n{config.upper()}:")
            report.append(f"  Average evaluations/sec: {summary.get('avg_evaluations_per_second', 0):.0f}")
            report.append(f"  Average time per eval: {summary.get('avg_time_per_evaluation_ns', 0):.1f} ns")
            report.append(f"  Scalability score: {summary.get('scalability_score', 0):.3f}")

        # Accuracy results
        if results.get('accuracy_results'):
            report.append("\\nAccuracy Analysis:")
            report.append("-" * 18)

            for config, accuracy in results['accuracy_results'].items():
                report.append(f"\\n{config.upper()}:")
                report.append(f"  Field at origin: {accuracy.get('field_at_origin_T', 0):.4f} T")
                report.append(f"  Cylindrical asymmetry: {accuracy.get('cylindrical_asymmetry', 0):.6f}")
                report.append(f"  Field gradient: {accuracy.get('field_gradient_T_per_cm', 0):.6f} T/cm")

        # Recommendations
        report.append("\\nRecommendations:")
        report.append("-" * 15)

        if results.get('performance_summary'):
            best_performance = max(results['performance_summary'].items(),
                                 key=lambda x: x[1].get('avg_evaluations_per_second', 0))
            report.append(f"• Best performance: {best_performance[0]} ({best_performance[1].get('avg_evaluations_per_second', 0):.0f} eval/s)")

            most_stable = min(results['performance_summary'].items(),
                            key=lambda x: x[1].get('scalability_score', float('inf')))
            report.append(f"• Most stable: {most_stable[0]} (scalability score: {most_stable[1].get('scalability_score', 0):.3f})")

        report_text = "\\n".join(report)

        # Save report
        report_file = self.output_dir / "benchmark_report.txt"
        with open(report_file, 'w') as f:
            f.write(report_text)

        logger.info(f"Report saved to {report_file}")
        return report_text


def main():
    parser = argparse.ArgumentParser(description='EPIC Field Performance Benchmark')
    parser.add_argument('--detector-path', default='/workspaces/epic',
                       help='Path to EPIC detector repository')
    parser.add_argument('--output-dir', default='benchmark_results',
                       help='Output directory for results')
    parser.add_argument('--config', choices=['marco_solenoid', 'lumi_magnets', 'all'],
                       default='all', help='Field configuration to test')
    parser.add_argument('--samples', type=int, nargs='+',
                       default=[1000, 10000, 100000],
                       help='Number of sample points to test')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Verify environment
    if 'DD4hepINSTALL' not in os.environ:
        logger.error("DD4hepINSTALL environment variable not set. Please source the EIC environment.")
        sys.exit(1)

    detector_path = Path(args.detector_path)
    if not detector_path.exists():
        logger.error(f"Detector path does not exist: {detector_path}")
        sys.exit(1)

    # Run benchmark
    benchmark = FieldBenchmark(detector_path, args.output_dir)
    benchmark.n_samples = args.samples

    if args.config != 'all':
        # Filter to specific configuration
        benchmark.field_configs = {args.config: benchmark.field_configs[args.config]}

    try:
        results = benchmark.run_comprehensive_benchmark()

        # Save results and generate report
        benchmark.save_results(results)
        report = benchmark.generate_report(results)

        print("\\nBenchmark Complete!")
        print("===================")
        print(report)

    except KeyboardInterrupt:
        logger.info("Benchmark interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Benchmark failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
