"""
Visualization and Reporting Module

This module provides comprehensive visualization and reporting capabilities including:
- Interactive plots using matplotlib/plotly
- HTML reports with summary statistics
- Sequence composition plots
- Quality score visualizations
- Variant analysis charts
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
import logging
from jinja2 import Template
import base64
import io


class Visualizer:
    """Main class for creating visualizations and reports."""
    
    def __init__(self, style: str = "seaborn-v0_8", figsize: Tuple[int, int] = (10, 6)):
        """
        Initialize the Visualizer.
        
        Args:
            style: Matplotlib style
            figsize: Default figure size
        """
        self.logger = logging.getLogger(__name__)
        plt.style.use(style)
        self.figsize = figsize
        sns.set_palette("husl")
    
    def plot_sequence_composition(self, composition_data: Dict[str, List[float]], 
                                output_file: Optional[str] = None, interactive: bool = False):
        """
        Plot sequence composition along the sequence.
        
        Args:
            composition_data: Dictionary with positions and composition data
            output_file: Output file path
            interactive: Whether to create interactive plot
        """
        if interactive:
            return self._plot_composition_interactive(composition_data, output_file)
        else:
            return self._plot_composition_static(composition_data, output_file)
    
    def _plot_composition_static(self, composition_data: Dict[str, List[float]], 
                               output_file: Optional[str] = None):
        """Create static composition plot with matplotlib."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=self.figsize, sharex=True)
        
        positions = composition_data["positions"]
        
        # GC content plot
        ax1.plot(positions, composition_data["gc_content"], 'b-', linewidth=2, label='GC content')
        ax1.fill_between(positions, composition_data["gc_content"], alpha=0.3, color='blue')
        ax1.set_ylabel('GC Content (%)')
        ax1.set_title('Sequence Composition Analysis')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # AT content plot
        ax2.plot(positions, composition_data["at_content"], 'r-', linewidth=2, label='AT content')
        ax2.fill_between(positions, composition_data["at_content"], alpha=0.3, color='red')
        ax2.set_ylabel('AT Content (%)')
        ax2.set_xlabel('Position (bp)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
    
    def _plot_composition_interactive(self, composition_data: Dict[str, List[float]], 
                                    output_file: Optional[str] = None):
        """Create interactive composition plot with plotly."""
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('GC Content', 'AT Content'),
            shared_xaxes=True,
            vertical_spacing=0.1
        )
        
        positions = composition_data["positions"]
        
        # GC content
        fig.add_trace(
            go.Scatter(
                x=positions,
                y=composition_data["gc_content"],
                mode='lines',
                name='GC Content',
                line=dict(color='blue', width=2),
                fill='tozeroy',
                fillcolor='rgba(0,0,255,0.3)'
            ),
            row=1, col=1
        )
        
        # AT content
        fig.add_trace(
            go.Scatter(
                x=positions,
                y=composition_data["at_content"],
                mode='lines',
                name='AT Content',
                line=dict(color='red', width=2),
                fill='tozeroy',
                fillcolor='rgba(255,0,0,0.3)'
            ),
            row=2, col=1
        )
        
        fig.update_layout(
            title='Sequence Composition Analysis',
            xaxis2_title='Position (bp)',
            yaxis_title='GC Content (%)',
            yaxis2_title='AT Content (%)',
            height=600
        )
        
        if output_file:
            fig.write_html(output_file)
        
        return fig
    
    def plot_quality_scores(self, quality_scores: List[int], 
                           output_file: Optional[str] = None):
        """
        Plot quality score distribution.
        
        Args:
            quality_scores: List of quality scores
            output_file: Output file path
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Quality score histogram
        ax1.hist(quality_scores, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        ax1.axvline(np.mean(quality_scores), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(quality_scores):.1f}')
        ax1.axvline(20, color='orange', linestyle='--', label='Q20 threshold')
        ax1.axvline(30, color='green', linestyle='--', label='Q30 threshold')
        ax1.set_xlabel('Quality Score')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Quality Score Distribution')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Quality score along sequence
        positions = range(len(quality_scores))
        ax2.plot(positions, quality_scores, alpha=0.6, linewidth=0.5)
        ax2.axhline(20, color='orange', linestyle='--', label='Q20')
        ax2.axhline(30, color='green', linestyle='--', label='Q30')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Quality Score')
        ax2.set_title('Quality Scores Along Sequence')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_orfs(self, orfs: List[Dict], sequence_length: int, 
                  output_file: Optional[str] = None):
        """
        Plot Open Reading Frames.
        
        Args:
            orfs: List of ORF dictionaries
            sequence_length: Length of the sequence
            output_file: Output file path
        """
        fig, ax = plt.subplots(figsize=(15, 8))
        
        # Group ORFs by strand and frame
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']
        
        for i, orf in enumerate(orfs[:20]):  # Show top 20 ORFs
            y_pos = orf['strand'] * 3 + orf['frame']
            color = colors[(orf['frame'] - 1) % len(colors)]
            
            # Draw ORF as rectangle
            rect = patches.Rectangle(
                (orf['start'], y_pos - 0.4), 
                orf['length'], 0.8,
                linewidth=1, edgecolor='black', 
                facecolor=color, alpha=0.7
            )
            ax.add_patch(rect)
            
            # Add ORF label
            ax.text(orf['start'] + orf['length']/2, y_pos, 
                   f"ORF{i+1}\n{orf['length']}bp", 
                   ha='center', va='center', fontsize=8)
        
        # Set up plot
        ax.set_xlim(0, sequence_length)
        ax.set_ylim(-3.5, 3.5)
        ax.set_xlabel('Position (bp)')
        ax.set_ylabel('Reading Frame')
        ax.set_title('Open Reading Frames')
        
        # Add frame labels
        ax.set_yticks([-3, -2, -1, 1, 2, 3])
        ax.set_yticklabels(['-3', '-2', '-1', '+1', '+2', '+3'])
        ax.axhline(0, color='black', linewidth=1)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_variant_statistics(self, variants, output_file: Optional[str] = None):
        """
        Plot variant analysis statistics.
        
        Args:
            variants: List of Variant objects
            output_file: Output file path
        """
        if not variants:
            self.logger.warning("No variants to plot")
            return None
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # Variant type distribution
        variant_types = []
        for var in variants:
            if var.is_snp():
                variant_types.append('SNP')
            elif var.is_insertion():
                variant_types.append('Insertion')
            elif var.is_deletion():
                variant_types.append('Deletion')
            else:
                variant_types.append('Other')
        
        type_counts = pd.Series(variant_types).value_counts()
        ax1.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%')
        ax1.set_title('Variant Type Distribution')
        
        # Quality score distribution
        qualities = [var.quality for var in variants]
        ax2.hist(qualities, bins=30, alpha=0.7, color='lightblue', edgecolor='black')
        ax2.axvline(np.mean(qualities), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(qualities):.1f}')
        ax2.set_xlabel('Quality Score')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Variant Quality Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Allele frequency distribution
        allele_freqs = [var.info.get('AF', 0) for var in variants]
        ax3.hist(allele_freqs, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
        ax3.set_xlabel('Allele Frequency')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Allele Frequency Distribution')
        ax3.grid(True, alpha=0.3)
        
        # Variant positions
        positions = [var.pos for var in variants]
        ax4.scatter(positions, qualities, alpha=0.6, s=30)
        ax4.set_xlabel('Position')
        ax4.set_ylabel('Quality Score')
        ax4.set_title('Variant Quality vs Position')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
        
        return fig
    
    def create_html_report(self, analysis_results: Dict, output_file: str):
        """
        Create comprehensive HTML report.
        
        Args:
            analysis_results: Dictionary containing all analysis results
            output_file: Output HTML file path
        """
        # HTML template
        html_template = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>GenomicsToolkit Analysis Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
                .section { margin: 30px 0; }
                .stats-table { border-collapse: collapse; width: 100%; }
                .stats-table th, .stats-table td { 
                    border: 1px solid #ddd; padding: 8px; text-align: left; 
                }
                .stats-table th { background-color: #f2f2f2; }
                .plot-container { text-align: center; margin: 20px 0; }
                .warning { color: #d9534f; }
                .success { color: #5cb85c; }
                .info { color: #5bc0de; }
            </style>
        </head>
        <body>
            <div class="header">
                <h1>GenomicsToolkit Analysis Report</h1>
                <p>Generated on: {{ timestamp }}</p>
                <p>Analysis type: {{ analysis_type }}</p>
            </div>
            
            {% if sequence_stats %}
            <div class="section">
                <h2>Sequence Statistics</h2>
                <table class="stats-table">
                    <tr><th>Metric</th><th>Value</th></tr>
                    {% for key, value in sequence_stats.items() %}
                    <tr><td>{{ key }}</td><td>{{ value }}</td></tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            {% if orf_data %}
            <div class="section">
                <h2>Open Reading Frames</h2>
                <p>Found {{ orf_data|length }} ORFs</p>
                <table class="stats-table">
                    <tr><th>ORF</th><th>Start</th><th>End</th><th>Length</th><th>Strand</th><th>Frame</th></tr>
                    {% for orf in orf_data[:10] %}
                    <tr>
                        <td>ORF {{ loop.index }}</td>
                        <td>{{ orf.start }}</td>
                        <td>{{ orf.end }}</td>
                        <td>{{ orf.length }}</td>
                        <td>{{ '+' if orf.strand == 1 else '-' }}</td>
                        <td>{{ orf.frame }}</td>
                    </tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            {% if variant_stats %}
            <div class="section">
                <h2>Variant Analysis</h2>
                <table class="stats-table">
                    <tr><th>Metric</th><th>Value</th></tr>
                    {% for key, value in variant_stats.items() %}
                    <tr><td>{{ key }}</td><td>{{ value }}</td></tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            {% if quality_stats %}
            <div class="section">
                <h2>Quality Control</h2>
                <table class="stats-table">
                    <tr><th>Metric</th><th>Value</th></tr>
                    {% for key, value in quality_stats.items() %}
                    <tr>
                        <td>{{ key }}</td>
                        <td class="{% if 'q30' in key|lower and value > 80 %}success{% elif 'q20' in key|lower and value > 90 %}success{% else %}info{% endif %}">
                            {{ value }}
                        </td>
                    </tr>
                    {% endfor %}
                </table>
            </div>
            {% endif %}
            
            <div class="section">
                <h2>Plots and Visualizations</h2>
                {% for plot_name, plot_data in plots.items() %}
                <div class="plot-container">
                    <h3>{{ plot_name }}</h3>
                    <img src="data:image/png;base64,{{ plot_data }}" style="max-width: 100%;">
                </div>
                {% endfor %}
            </div>
            
            <div class="section">
                <h2>Analysis Summary</h2>
                <ul>
                    {% if sequence_stats %}
                    <li>Analyzed sequence of {{ sequence_stats.get('length', 'unknown') }} bp</li>
                    <li>GC content: {{ "%.1f"|format(sequence_stats.get('gc_content', 0)) }}%</li>
                    {% endif %}
                    {% if orf_data %}
                    <li>Identified {{ orf_data|length }} open reading frames</li>
                    {% endif %}
                    {% if variant_stats %}
                    <li>Detected {{ variant_stats.get('total_variants', 0) }} variants</li>
                    {% endif %}
                </ul>
            </div>
            
            <div class="section">
                <p><em>Report generated by GenomicsToolkit v1.0.0</em></p>
            </div>
        </body>
        </html>
        """
        
        # Convert plots to base64 for embedding
        plots_b64 = {}
        if 'plots' in analysis_results:
            for plot_name, fig in analysis_results['plots'].items():
                img_buffer = io.BytesIO()
                fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                img_buffer.seek(0)
                plots_b64[plot_name] = base64.b64encode(img_buffer.getvalue()).decode()
                plt.close(fig)
        
        # Render template
        template = Template(html_template)
        html_content = template.render(
            timestamp=pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
            analysis_type=analysis_results.get('analysis_type', 'Complete Analysis'),
            sequence_stats=analysis_results.get('sequence_stats'),
            orf_data=analysis_results.get('orf_data'),
            variant_stats=analysis_results.get('variant_stats'),
            quality_stats=analysis_results.get('quality_stats'),
            plots=plots_b64
        )
        
        # Write HTML file
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"HTML report saved to {output_file}")
    
    def export_results(self, results: Dict, output_dir: str, formats: List[str] = ['csv', 'json']):
        """
        Export analysis results in various formats.
        
        Args:
            results: Analysis results dictionary
            output_dir: Output directory
            formats: List of export formats ('csv', 'json', 'vcf', 'bed')
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        for fmt in formats:
            if fmt == 'csv':
                self._export_csv(results, output_path)
            elif fmt == 'json':
                self._export_json(results, output_path)
            elif fmt == 'vcf' and 'variants' in results:
                self._export_vcf(results['variants'], output_path)
            elif fmt == 'bed' and 'orfs' in results:
                self._export_bed(results['orfs'], output_path)
    
    def _export_csv(self, results: Dict, output_path: Path):
        """Export results to CSV format."""
        if 'variants' in results:
            variant_data = []
            for var in results['variants']:
                variant_data.append({
                    'chromosome': var.chrom,
                    'position': var.pos,
                    'reference': var.ref,
                    'alternative': var.alt,
                    'quality': var.quality,
                    'type': 'SNP' if var.is_snp() else 'INDEL'
                })
            
            df = pd.DataFrame(variant_data)
            df.to_csv(output_path / 'variants.csv', index=False)
        
        if 'orfs' in results:
            orf_df = pd.DataFrame(results['orfs'])
            orf_df.to_csv(output_path / 'orfs.csv', index=False)
    
    def _export_json(self, results: Dict, output_path: Path):
        """Export results to JSON format."""
        import json
        
        # Convert results to JSON-serializable format
        json_results = {}
        for key, value in results.items():
            if key == 'variants':
                json_results[key] = [
                    {
                        'chrom': var.chrom,
                        'pos': var.pos,
                        'ref': var.ref,
                        'alt': var.alt,
                        'quality': var.quality,
                        'info': var.info
                    } for var in value
                ]
            else:
                json_results[key] = value
        
        with open(output_path / 'results.json', 'w') as f:
            json.dump(json_results, f, indent=2)
    
    def _export_vcf(self, variants, output_path: Path):
        """Export variants to VCF format."""
        from .variant_calling import VariantCaller
        caller = VariantCaller()
        caller.write_vcf(variants, str(output_path / 'variants.vcf'))
    
    def _export_bed(self, orfs: List[Dict], output_path: Path):
        """Export ORFs to BED format."""
        with open(output_path / 'orfs.bed', 'w') as f:
            for i, orf in enumerate(orfs):
                # BED format: chrom, start, end, name, score, strand
                strand = '+' if orf['strand'] == 1 else '-'
                f.write(f"chr1\t{orf['start']}\t{orf['end']}\tORF_{i+1}\t{orf['length']}\t{strand}\n")