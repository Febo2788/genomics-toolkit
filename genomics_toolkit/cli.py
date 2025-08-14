"""
Command Line Interface Module

This module provides a comprehensive CLI for the GenomicsToolkit using Click.
Supports all major operations including sequence analysis, variant calling, and visualization.
"""

import click
import logging
import sys
from pathlib import Path
from typing import Optional, List
import json

from .sequence_analysis import SequenceAnalyzer
from .variant_calling import VariantCaller
from .visualization import Visualizer
from .utils import FileHandler, QualityControl, ConfigManager, timing_decorator


# Set up logging
def setup_logging(log_level: str, log_file: Optional[str] = None):
    """Set up logging configuration."""
    level = getattr(logging, log_level.upper())
    
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


@click.group()
@click.option('--log-level', default='INFO', 
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='Set logging level')
@click.option('--log-file', help='Log file path')
@click.option('--config', help='Configuration file path')
@click.pass_context
def cli(ctx, log_level, log_file, config):
    """GenomicsToolkit: Comprehensive bioinformatics pipeline for sequence analysis and variant calling."""
    setup_logging(log_level, log_file)
    
    # Initialize context
    ctx.ensure_object(dict)
    ctx.obj['config'] = ConfigManager(config)
    ctx.obj['file_handler'] = FileHandler()
    ctx.obj['qc'] = QualityControl()


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', default='./output', help='Output directory')
@click.option('--format', 'file_format', help='Force file format (auto-detect if not specified)')
@click.option('--min-orf-length', default=100, help='Minimum ORF length in nucleotides')
@click.option('--window-size', default=100, help='Window size for composition analysis')
@click.option('--report', is_flag=True, help='Generate HTML report')
@click.option('--plots', is_flag=True, help='Generate plots')
@click.pass_context
@timing_decorator
def analyze_sequence(ctx, input_file, output_dir, file_format, min_orf_length, window_size, report, plots):
    """Analyze DNA/RNA sequences from FASTA or FASTQ files."""
    click.echo(f"Analyzing sequence file: {input_file}")
    
    # Initialize components
    analyzer = SequenceAnalyzer()
    visualizer = Visualizer()
    file_handler = ctx.obj['file_handler']
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Detect file format
    if not file_format:
        file_format = file_handler.detect_file_format(input_file)
    
    click.echo(f"Detected file format: {file_format}")
    
    results = {}
    
    try:
        if file_format == 'fasta':
            sequences = list(analyzer.parse_fasta(input_file))
            
            # Analyze each sequence
            all_stats = []
            all_orfs = []
            all_composition = []
            
            for i, (header, sequence) in enumerate(sequences):
                click.echo(f"Processing sequence {i+1}: {header}")
                
                # Basic statistics
                stats = analyzer.basic_stats(sequence)
                stats['header'] = header
                all_stats.append(stats)
                
                # ORF detection
                orfs = analyzer.find_orfs(sequence, min_orf_length)
                all_orfs.extend(orfs)
                
                # Composition analysis
                composition = analyzer.composition_analysis(sequence, window_size)
                all_composition.append(composition)
                
                # Generate plots if requested
                if plots and i == 0:  # Plot first sequence
                    comp_plot = visualizer.plot_sequence_composition(
                        composition, 
                        output_path / f"composition_{i+1}.png"
                    )
                    
                    if orfs:
                        orf_plot = visualizer.plot_orfs(
                            orfs, len(sequence),
                            output_path / f"orfs_{i+1}.png"
                        )
            
            results.update({
                'analysis_type': 'Sequence Analysis',
                'sequence_stats': all_stats,
                'orf_data': all_orfs,
                'composition_data': all_composition
            })
            
        elif file_format == 'fastq':
            sequences = list(analyzer.parse_fastq(input_file))
            
            # Analyze quality scores
            all_quality_stats = []
            all_sequences = []
            
            for i, (header, sequence, quality_scores) in enumerate(sequences):
                if i >= 1000:  # Limit for demo
                    break
                    
                # Quality analysis
                quality_stats = analyzer.analyze_quality_scores(quality_scores)
                quality_stats['header'] = header
                all_quality_stats.append(quality_stats)
                
                # Sequence statistics
                seq_stats = analyzer.basic_stats(sequence)
                seq_stats['header'] = header
                all_sequences.append(seq_stats)
                
                # Plot quality scores for first sequence
                if plots and i == 0:
                    quality_plot = visualizer.plot_quality_scores(
                        quality_scores,
                        output_path / f"quality_{i+1}.png"
                    )
            
            results.update({
                'analysis_type': 'FASTQ Quality Analysis',
                'sequence_stats': all_sequences,
                'quality_stats': all_quality_stats
            })
        
        # Export results
        if results:
            # Save JSON results
            with open(output_path / 'analysis_results.json', 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            # Generate HTML report if requested
            if report:
                visualizer.create_html_report(results, output_path / 'report.html')
                click.echo(f"HTML report generated: {output_path / 'report.html'}")
            
            # Summary statistics
            if 'sequence_stats' in results:
                total_sequences = len(results['sequence_stats'])
                total_length = sum(s['length'] for s in results['sequence_stats'])
                avg_gc = sum(s['gc_content'] for s in results['sequence_stats']) / total_sequences
                
                click.echo(f"\nAnalysis Summary:")
                click.echo(f"Total sequences: {total_sequences}")
                click.echo(f"Total length: {total_length:,} bp")
                click.echo(f"Average GC content: {avg_gc:.1f}%")
                
                if 'orf_data' in results:
                    click.echo(f"ORFs found: {len(results['orf_data'])}")
        
        click.echo(f"\nAnalysis complete! Results saved to {output_dir}")
        
    except Exception as e:
        click.echo(f"Error during analysis: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument('reference_file', type=click.Path(exists=True))
@click.option('--bam-file', type=click.Path(exists=True), help='BAM alignment file')
@click.option('--simulate', is_flag=True, help='Simulate variants for testing')
@click.option('--num-snps', default=100, help='Number of SNPs to simulate')
@click.option('--num-indels', default=20, help='Number of indels to simulate')
@click.option('--output-dir', '-o', default='./output', help='Output directory')
@click.option('--min-coverage', default=10, help='Minimum coverage depth')
@click.option('--min-quality', default=20.0, help='Minimum variant quality')
@click.option('--min-allele-freq', default=0.2, help='Minimum allele frequency')
@click.option('--report', is_flag=True, help='Generate HTML report')
@click.pass_context
@timing_decorator
def call_variants(ctx, reference_file, bam_file, simulate, num_snps, num_indels, 
                 output_dir, min_coverage, min_quality, min_allele_freq, report):
    """Call variants from alignment data or simulate variants for testing."""
    click.echo("Starting variant calling...")
    
    # Initialize components
    caller = VariantCaller()
    visualizer = Visualizer()
    analyzer = SequenceAnalyzer()
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    try:
        if simulate:
            click.echo("Simulating variants for testing...")
            
            # Read reference sequence
            sequences = list(analyzer.parse_fasta(reference_file))
            if not sequences:
                click.echo("No sequences found in reference file", err=True)
                sys.exit(1)
            
            ref_sequence = sequences[0][1]  # Use first sequence
            
            # Simulate variants
            variants = caller.simulate_variants(ref_sequence, num_snps, num_indels)
            
        elif bam_file:
            click.echo(f"Calling variants from BAM file: {bam_file}")
            variants = caller.call_variants_from_alignment(
                bam_file, reference_file, min_coverage, min_quality, min_allele_freq
            )
        else:
            click.echo("Either --bam-file or --simulate must be specified", err=True)
            sys.exit(1)
        
        if not variants:
            click.echo("No variants found", err=True)
            return
        
        click.echo(f"Found {len(variants)} variants")
        
        # Filter variants
        filtered_variants = caller.filter_variants(
            min_quality=min_quality,
            min_allele_freq=min_allele_freq
        )
        
        click.echo(f"After filtering: {len(filtered_variants)} variants")
        
        # Annotate variants
        annotated_variants = caller.annotate_variants(filtered_variants)
        
        # Calculate population statistics
        pop_stats = caller.calculate_population_stats(annotated_variants)
        
        # Write VCF file
        caller.write_vcf(annotated_variants, output_path / 'variants.vcf')
        
        # Generate plots
        variant_plot = visualizer.plot_variant_statistics(
            annotated_variants, 
            output_path / 'variant_statistics.png'
        )
        
        # Prepare results for report
        results = {
            'analysis_type': 'Variant Calling',
            'variants': annotated_variants,
            'variant_stats': pop_stats,
            'plots': {'Variant Statistics': variant_plot}
        }
        
        # Generate HTML report if requested
        if report:
            visualizer.create_html_report(results, output_path / 'variant_report.html')
            click.echo(f"HTML report generated: {output_path / 'variant_report.html'}")
        
        # Print summary
        click.echo(f"\nVariant Calling Summary:")
        click.echo(f"Total variants: {pop_stats['total_variants']}")
        click.echo(f"SNPs: {pop_stats['snp_count']}")
        click.echo(f"Indels: {pop_stats['indel_count']}")
        click.echo(f"Ti/Tv ratio: {pop_stats['ti_tv_ratio']:.2f}")
        
        click.echo(f"\nResults saved to {output_dir}")
        
    except Exception as e:
        click.echo(f"Error during variant calling: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', default='./output', help='Output directory')
@click.option('--interactive', is_flag=True, help='Generate interactive plots')
@click.option('--format', 'plot_format', default='png', 
              type=click.Choice(['png', 'pdf', 'svg', 'html']),
              help='Output format for plots')
@click.pass_context
def visualize(ctx, input_file, output_dir, interactive, plot_format):
    """Generate visualizations from analysis results."""
    click.echo(f"Creating visualizations from: {input_file}")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    visualizer = Visualizer()
    
    try:
        # Load results
        with open(input_file, 'r') as f:
            results = json.load(f)
        
        # Generate appropriate visualizations based on data type
        if 'composition_data' in results:
            for i, comp_data in enumerate(results['composition_data']):
                visualizer.plot_sequence_composition(
                    comp_data,
                    output_path / f"composition_{i+1}.{plot_format}",
                    interactive=interactive
                )
        
        if 'quality_stats' in results and results['quality_stats']:
            # Create mock quality scores for visualization
            quality_scores = [30] * 100  # Placeholder
            visualizer.plot_quality_scores(
                quality_scores,
                output_path / f"quality_scores.{plot_format}"
            )
        
        if 'variants' in results:
            # Convert variant dictionaries back to Variant objects
            from .variant_calling import Variant
            variants = []
            for v_data in results['variants']:
                variant = Variant(v_data['chrom'], v_data['pos'], v_data['ref'], v_data['alt'])
                variant.quality = v_data['quality']
                variant.info = v_data['info']
                variants.append(variant)
            
            visualizer.plot_variant_statistics(
                variants,
                output_path / f"variants.{plot_format}"
            )
        
        click.echo(f"Visualizations saved to {output_dir}")
        
    except Exception as e:
        click.echo(f"Error creating visualizations: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--key', help='Configuration key to get/set')
@click.option('--value', help='Value to set (if not provided, will get the value)')
@click.option('--list-all', is_flag=True, help='List all configuration settings')
@click.pass_context
def config(ctx, key, value, list_all):
    """Manage configuration settings."""
    config_manager = ctx.obj['config']
    
    if list_all:
        click.echo("Current configuration:")
        click.echo(json.dumps(config_manager.config, indent=2))
    elif key and value:
        config_manager.set(key, value)
        click.echo(f"Set {key} = {value}")
    elif key:
        current_value = config_manager.get(key)
        click.echo(f"{key} = {current_value}")
    else:
        click.echo("Use --key to get a value, --key and --value to set, or --list-all to see all settings")


@cli.command()
@click.argument('file_path', type=click.Path(exists=True))
@click.pass_context
def validate(ctx, file_path):
    """Validate input files and check system requirements."""
    file_handler = ctx.obj['file_handler']
    qc = ctx.obj['qc']
    
    click.echo(f"Validating file: {file_path}")
    
    # File validation
    validation_result = file_handler.validate_file_integrity(file_path)
    
    if validation_result['valid']:
        click.echo("✓ File is valid")
        click.echo(f"  Format: {validation_result['file_format']}")
        click.echo(f"  Size: {validation_result['file_size']:,} bytes")
        click.echo(f"  Lines: {validation_result['line_count']:,}")
    else:
        click.echo(f"✗ File validation failed: {validation_result['error']}")
        return
    
    # System resource check
    resources = qc.check_system_resources()
    click.echo(f"\nSystem Resources:")
    click.echo(f"  CPU cores: {resources['cpu_count']}")
    click.echo(f"  Available memory: {resources['memory_available'] / 1e9:.1f} GB")
    click.echo(f"  Free disk space: {resources['disk_free'] / 1e9:.1f} GB")
    
    # Estimate processing time
    processing_time = qc.estimate_processing_time(validation_result['file_size'])
    click.echo(f"  Estimated processing time: {processing_time:.1f} seconds")


@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', default='./pipeline_output', help='Output directory')
@click.option('--reference', type=click.Path(exists=True), help='Reference genome file')
@click.option('--skip-variants', is_flag=True, help='Skip variant calling')
@click.option('--threads', default=4, help='Number of threads to use')
@click.pass_context
@timing_decorator
def pipeline(ctx, input_file, output_dir, reference, skip_variants, threads):
    """Run complete genomics analysis pipeline."""
    click.echo("Starting complete genomics analysis pipeline...")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    file_handler = ctx.obj['file_handler']
    file_format = file_handler.detect_file_format(input_file)
    
    try:
        # Step 1: Sequence Analysis
        click.echo("Step 1: Sequence analysis...")
        ctx.invoke(analyze_sequence, 
                  input_file=input_file,
                  output_dir=str(output_path / 'sequence_analysis'),
                  plots=True,
                  report=True)
        
        # Step 2: Variant Calling (if reference provided and not skipped)
        if reference and not skip_variants:
            click.echo("Step 2: Variant calling...")
            ctx.invoke(call_variants,
                      reference_file=reference,
                      simulate=True,  # Use simulation if no BAM file
                      output_dir=str(output_path / 'variant_calling'),
                      report=True)
        
        # Step 3: Generate comprehensive report
        click.echo("Step 3: Generating comprehensive report...")
        
        # Combine all results
        combined_results = {
            'analysis_type': 'Complete Pipeline',
            'input_file': str(input_file),
            'file_format': file_format,
            'timestamp': str(pd.Timestamp.now())
        }
        
        # Save pipeline summary
        with open(output_path / 'pipeline_summary.json', 'w') as f:
            json.dump(combined_results, f, indent=2)
        
        click.echo(f"\n🎉 Pipeline completed successfully!")
        click.echo(f"All results saved to: {output_dir}")
        click.echo(f"Check the following directories:")
        click.echo(f"  - Sequence analysis: {output_path / 'sequence_analysis'}")
        if reference and not skip_variants:
            click.echo(f"  - Variant calling: {output_path / 'variant_calling'}")
        
    except Exception as e:
        click.echo(f"Pipeline failed: {e}", err=True)
        sys.exit(1)


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == '__main__':
    main()