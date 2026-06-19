#!/usr/bin/env python3
"""
CheckM Quality Assessment for NXTTHRAX Pipeline

Runs CheckM on all input genomes, generates a contamination/completeness report,
filters genomes by quality thresholds, and optionally updates metadata.tsv.

Usage:
    python run_checkm_qc.py \
        --genome-dir data/genomes \
        --output-dir checkm_results \
        --metadata config/metadata.tsv \
        --min-completeness 95 \
        --max-contamination 5 \
        --update-metadata
"""

import argparse
import subprocess
import sys
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


def check_checkm_installed():
    """Verify CheckM is installed and accessible."""
    try:
        result = subprocess.run(['checkm', '--version'], capture_output=True, text=True)
        logger.info(f"CheckM found: {result.stdout.strip()}")
        return True
    except FileNotFoundError:
        logger.error("CheckM not found. Install with: conda install -c bioconda checkm-genome")
        return False


def get_genome_files(genome_dir):
    """Find all genome files in directory."""
    genome_dir = Path(genome_dir)
    extensions = ['*.fna', '*.fasta', '*.fa', '*.fna.gz', '*.fasta.gz', '*.fa.gz']
    
    genomes = []
    for ext in extensions:
        genomes.extend(genome_dir.glob(ext))
    
    if not genomes:
        logger.error(f"No genome files found in {genome_dir}")
        return []
    
    logger.info(f"Found {len(genomes)} genome file(s)")
    return sorted(genomes)


def run_checkm(genome_dir, output_dir, threads=4):
    """Run CheckM on all genomes."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Running CheckM on all genomes...")
    logger.info(f"Output directory: {output_dir}")
    
    try:
        cmd = [
            'checkm', 'lineage_wf',
            '-x', 'fna',  # or fasta/fa depending on your files
            '-t', str(threads),
            '--json',
            str(genome_dir),
            str(output_dir)
        ]
        
        logger.info(f"Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"CheckM failed with return code {result.returncode}")
            logger.error(f"stderr: {result.stderr}")
            return False
        
        logger.info("CheckM completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error running CheckM: {e}")
        return False


def parse_checkm_output(output_dir):
    """Parse CheckM output and extract metrics."""
    output_dir = Path(output_dir)
    results_file = output_dir / 'results.json'
    
    if not results_file.exists():
        logger.error(f"CheckM results file not found: {results_file}")
        return None
    
    try:
        with open(results_file, 'r') as f:
            results = json.load(f)
        
        logger.info(f"Parsed results for {len(results)} genomes")
        return results
        
    except Exception as e:
        logger.error(f"Error parsing CheckM output: {e}")
        return None


def generate_report(results, output_file='checkm_quality_report.tsv'):
    """Generate a human-readable quality report."""
    report_data = []
    
    for genome_id, data in results.items():
        report_data.append({
            'genome': genome_id,
            'completeness': data.get('Completeness', 'N/A'),
            'contamination': data.get('Contamination', 'N/A'),
            'strain_heterogeneity': data.get('Strain heterogeneity', 'N/A'),
            'genome_size': data.get('Genome size (bp)', 'N/A'),
            'gc': data.get('GC (%)', 'N/A'),
            'num_contigs': data.get('# contigs', 'N/A'),
        })
    
    df = pd.DataFrame(report_data)
    df = df.sort_values('contamination').reset_index(drop=True)
    
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Report written to: {output_file}")
    
    return df


def filter_genomes(results, min_completeness=95, max_contamination=5):
    """Filter genomes by quality thresholds."""
    passed = {}
    failed = {}
    
    for genome_id, data in results.items():
        completeness = float(data.get('Completeness', 0))
        contamination = float(data.get('Contamination', 100))
        
        if completeness >= min_completeness and contamination <= max_contamination:
            passed[genome_id] = data
        else:
            failed[genome_id] = {
                'completeness': completeness,
                'contamination': contamination,
                'reason': []
            }
            if completeness < min_completeness:
                failed[genome_id]['reason'].append(
                    f"low completeness ({completeness:.1f}% < {min_completeness}%)"
                )
            if contamination > max_contamination:
                failed[genome_id]['reason'].append(
                    f"high contamination ({contamination:.1f}% > {max_contamination}%)"
                )
    
    logger.info(f"Quality filtering results:")
    logger.info(f"  Passed: {len(passed)} genomes")
    logger.info(f"  Failed: {len(failed)} genomes")
    
    return passed, failed


def write_failed_genomes(failed_genomes, output_file='genomes_failed_qc.tsv'):
    """Write list of failed genomes with reasons."""
    failed_data = []
    
    for genome_id, data in failed_genomes.items():
        failed_data.append({
            'genome': genome_id,
            'completeness': data['completeness'],
            'contamination': data['contamination'],
            'reason': '; '.join(data['reason'])
        })
    
    if failed_data:
        df = pd.DataFrame(failed_data)
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Failed genomes written to: {output_file}")
    else:
        logger.info("All genomes passed quality filtering!")


def update_metadata(metadata_file, results, output_file=None):
    """Update metadata.tsv with CheckM metrics."""
    if not Path(metadata_file).exists():
        logger.warning(f"Metadata file not found: {metadata_file}")
        return
    
    if output_file is None:
        output_file = metadata_file
    
    try:
        # Read existing metadata
        df = pd.read_csv(metadata_file, sep='\t')
        
        # Create mapping from genome filename to CheckM metrics
        checkm_dict = {}
        for genome_id, data in results.items():
            # genome_id might be genome name without extension
            # Try to match with strain column
            checkm_dict[genome_id] = {
                'checkm_completeness': round(float(data.get('Completeness', 0)), 2),
                'checkm_contamination': round(float(data.get('Contamination', 0)), 2),
                'checkm_strain_heterogeneity': round(float(data.get('Strain heterogeneity', 0)), 2),
            }
        
        # Update metadata
        updated_rows = 0
        for idx, row in df.iterrows():
            strain = row['strain']
            # Remove common extensions
            strain_base = strain.replace('.fna', '').replace('.fasta', '').replace('.fa', '')
            
            if strain_base in checkm_dict or strain in checkm_dict:
                key = strain_base if strain_base in checkm_dict else strain
                metrics = checkm_dict[key]
                
                # Update or create columns
                for col, value in metrics.items():
                    df.loc[idx, col] = value
                
                updated_rows += 1
        
        # Write updated metadata
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Updated metadata written to: {output_file}")
        logger.info(f"Updated {updated_rows} rows with CheckM metrics")
        
    except Exception as e:
        logger.error(f"Error updating metadata: {e}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--genome-dir',
        required=True,
        help='Directory containing genome files (.fna, .fasta, .fa)'
    )
    
    parser.add_argument(
        '--output-dir',
        default='checkm_results',
        help='Directory for CheckM output (default: checkm_results)'
    )
    
    parser.add_argument(
        '--metadata',
        default=None,
        help='Path to metadata.tsv file to update with CheckM metrics'
    )
    
    parser.add_argument(
        '--min-completeness',
        type=float,
        default=95,
        help='Minimum completeness threshold (default: 95%%)'
    )
    
    parser.add_argument(
        '--max-contamination',
        type=float,
        default=5,
        help='Maximum contamination threshold (default: 5%%)'
    )
    
    parser.add_argument(
        '--update-metadata',
        action='store_true',
        help='Update metadata.tsv with CheckM metrics'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for CheckM (default: 4)'
    )
    
    args = parser.parse_args()
    
    # Verify CheckM is installed
    if not check_checkm_installed():
        sys.exit(1)
    
    # Check genome files exist
    genomes = get_genome_files(args.genome_dir)
    if not genomes:
        sys.exit(1)
    
    # Run CheckM
    if not run_checkm(args.genome_dir, args.output_dir, threads=args.threads):
        sys.exit(1)
    
    # Parse results
    results = parse_checkm_output(args.output_dir)
    if results is None:
        sys.exit(1)
    
    # Generate report
    report_df = generate_report(results)
    logger.info("\n" + "="*60)
    logger.info("CheckM Quality Summary:")
    logger.info("="*60)
    logger.info(f"\nCompleteness (mean): {report_df['completeness'].mean():.2f}%")
    logger.info(f"Contamination (mean): {report_df['contamination'].mean():.2f}%")
    logger.info(f"\nCompleteness range: {report_df['completeness'].min():.2f}% - {report_df['completeness'].max():.2f}%")
    logger.info(f"Contamination range: {report_df['contamination'].min():.2f}% - {report_df['contamination'].max():.2f}%")
    logger.info("="*60 + "\n")
    
    # Filter genomes
    passed, failed = filter_genomes(
        results,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination
    )
    write_failed_genomes(failed)
    
    # Update metadata if requested
    if args.update_metadata and args.metadata:
        update_metadata(args.metadata, results)
    
    logger.info("\n✓ CheckM quality assessment complete!")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Report file: checkm_quality_report.tsv")
    if failed:
        logger.info(f"Failed genomes: genomes_failed_qc.tsv")


if __name__ == '__main__':
    main()
