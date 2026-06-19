#!/usr/bin/env python3
"""
BUSCO Quality Assessment for NXTTHRAX Pipeline

Runs BUSCO on all input genomes, generates a completeness report,
filters genomes by quality thresholds, and optionally updates metadata.tsv.

BUSCO (Benchmarking Universal Single-Copy Orthologs) assesses genome 
completeness by searching for conserved single-copy orthologs.

Usage:
    python run_busco_qc.py \
        --genome-dir data/genomes \
        --output-dir busco_results \
        --metadata config/metadata.tsv \
        --lineage bacillales_odb10 \
        --min-completeness 90 \
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
import tempfile
import shutil

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


def check_busco_installed():
    """Verify BUSCO is installed and accessible."""
    try:
        result = subprocess.run(['busco', '--version'], capture_output=True, text=True)
        logger.info(f"BUSCO found: {result.stdout.strip()}")
        return True
    except FileNotFoundError:
        logger.error("BUSCO not found. Install with: conda install -c bioconda busco")
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


def run_busco_on_genome(genome_path, output_dir, lineage='bacillales_odb10', threads=4, mode='genome'):
    """Run BUSCO on a single genome."""
    genome_name = genome_path.stem.replace('.fna', '').replace('.fasta', '').replace('.fa', '')
    
    try:
        cmd = [
            'busco',
            '-i', str(genome_path),
            '-l', lineage,
            '-o', genome_name,
            '-m', mode,
            '-c', str(threads),
            '--out_path', str(output_dir),
            '--quiet'
        ]
        
        logger.debug(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        
        if result.returncode != 0:
            logger.warning(f"BUSCO failed for {genome_name}: {result.stderr}")
            return None
        
        return genome_name
        
    except subprocess.TimeoutExpired:
        logger.warning(f"BUSCO timeout for {genome_name} (>1 hour)")
        return None
    except Exception as e:
        logger.warning(f"Error running BUSCO on {genome_name}: {e}")
        return None


def parse_busco_results(output_dir, lineage='bacillales_odb10'):
    """Parse BUSCO output and extract metrics."""
    output_dir = Path(output_dir)
    results = {}
    
    # Find all genome result directories
    for genome_dir in sorted(output_dir.iterdir()):
        if not genome_dir.is_dir():
            continue
        
        genome_name = genome_dir.name
        summary_file = genome_dir / 'short_summary.json'
        
        if not summary_file.exists():
            logger.debug(f"Summary file not found for {genome_name}")
            continue
        
        try:
            with open(summary_file, 'r') as f:
                summary = json.load(f)
            
            # Extract metrics from BUSCO summary
            results_dict = summary.get('results', {})
            
            results[genome_name] = {
                'complete_single': results_dict.get('single_copy_number', 0),
                'complete_duplicated': results_dict.get('duplicated_number', 0),
                'fragmented': results_dict.get('fragmented_number', 0),
                'missing': results_dict.get('missing_number', 0),
                'total_length': results_dict.get('total_length', 0),
                'completeness': results_dict.get('complete_percent', 0.0),
            }
            
        except Exception as e:
            logger.warning(f"Error parsing BUSCO results for {genome_name}: {e}")
    
    if results:
        logger.info(f"Parsed results for {len(results)} genomes")
    else:
        logger.warning("No BUSCO results found")
    
    return results


def generate_report(results, output_file='busco_quality_report.tsv'):
    """Generate a human-readable BUSCO report."""
    report_data = []
    
    for genome_id, data in results.items():
        report_data.append({
            'genome': genome_id,
            'completeness': round(data.get('completeness', 0.0), 2),
            'complete_single_copy': data.get('complete_single', 0),
            'complete_duplicated': data.get('complete_duplicated', 0),
            'fragmented': data.get('fragmented', 0),
            'missing': data.get('missing', 0),
        })
    
    df = pd.DataFrame(report_data)
    df = df.sort_values('completeness', ascending=False).reset_index(drop=True)
    
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Report written to: {output_file}")
    
    return df


def filter_genomes(results, min_completeness=90):
    """Filter genomes by quality thresholds."""
    passed = {}
    failed = {}
    
    for genome_id, data in results.items():
        completeness = float(data.get('completeness', 0))
        
        if completeness >= min_completeness:
            passed[genome_id] = data
        else:
            failed[genome_id] = {
                'completeness': completeness,
                'reason': f"low completeness ({completeness:.1f}% < {min_completeness}%)"
            }
    
    logger.info(f"Quality filtering results:")
    logger.info(f"  Passed: {len(passed)} genomes")
    logger.info(f"  Failed: {len(failed)} genomes")
    
    return passed, failed


def write_failed_genomes(failed_genomes, output_file='genomes_failed_busco.tsv'):
    """Write list of failed genomes with reasons."""
    failed_data = []
    
    for genome_id, data in failed_genomes.items():
        failed_data.append({
            'genome': genome_id,
            'completeness': data['completeness'],
            'reason': data['reason']
        })
    
    if failed_data:
        df = pd.DataFrame(failed_data)
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Failed genomes written to: {output_file}")
    else:
        logger.info("All genomes passed BUSCO quality filtering!")


def update_metadata(metadata_file, results, output_file=None):
    """Update metadata.tsv with BUSCO metrics."""
    if not Path(metadata_file).exists():
        logger.warning(f"Metadata file not found: {metadata_file}")
        return
    
    if output_file is None:
        output_file = metadata_file
    
    try:
        # Read existing metadata
        df = pd.read_csv(metadata_file, sep='\t')
        
        # Create mapping from genome filename to BUSCO metrics
        busco_dict = {}
        for genome_id, data in results.items():
            busco_dict[genome_id] = {
                'busco_completeness': round(float(data.get('completeness', 0)), 2),
                'busco_complete_single': int(data.get('complete_single', 0)),
                'busco_complete_duplicated': int(data.get('complete_duplicated', 0)),
                'busco_fragmented': int(data.get('fragmented', 0)),
                'busco_missing': int(data.get('missing', 0)),
            }
        
        # Update metadata
        updated_rows = 0
        for idx, row in df.iterrows():
            strain = row['strain']
            # Remove common extensions
            strain_base = strain.replace('.fna', '').replace('.fasta', '').replace('.fa', '')
            
            if strain_base in busco_dict or strain in busco_dict:
                key = strain_base if strain_base in busco_dict else strain
                metrics = busco_dict[key]
                
                # Update or create columns
                for col, value in metrics.items():
                    df.loc[idx, col] = value
                
                updated_rows += 1
        
        # Write updated metadata
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Updated metadata written to: {output_file}")
        logger.info(f"Updated {updated_rows} rows with BUSCO metrics")
        
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
        default='busco_results',
        help='Directory for BUSCO output (default: busco_results)'
    )
    
    parser.add_argument(
        '--metadata',
        default=None,
        help='Path to metadata.tsv file to update with BUSCO metrics'
    )
    
    parser.add_argument(
        '--lineage',
        default='bacillales_odb10',
        help='BUSCO lineage dataset (default: bacillales_odb10)'
    )
    
    parser.add_argument(
        '--min-completeness',
        type=float,
        default=90,
        help='Minimum completeness threshold in %% (default: 90)'
    )
    
    parser.add_argument(
        '--update-metadata',
        action='store_true',
        help='Update metadata.tsv with BUSCO metrics'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads for BUSCO (default: 4)'
    )
    
    args = parser.parse_args()
    
    # Verify BUSCO is installed
    if not check_busco_installed():
        sys.exit(1)
    
    # Check genome files exist
    genomes = get_genome_files(args.genome_dir)
    if not genomes:
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run BUSCO on each genome
    logger.info(f"Running BUSCO on {len(genomes)} genome(s) with lineage: {args.lineage}")
    logger.info(f"(This may take 10-60 minutes depending on genome size and CPU cores)")
    logger.info(f"Output directory: {output_dir}")
    
    successful = 0
    for i, genome_path in enumerate(genomes, 1):
        logger.info(f"[{i}/{len(genomes)}] Processing {genome_path.name}...")
        result = run_busco_on_genome(genome_path, output_dir, args.lineage, args.threads)
        if result:
            successful += 1
    
    logger.info(f"BUSCO completed successfully on {successful}/{len(genomes)} genomes")
    
    # Parse results
    results = parse_busco_results(output_dir, args.lineage)
    if not results:
        logger.warning("No BUSCO results to process")
        sys.exit(1)
    
    # Generate report
    report_df = generate_report(results)
    logger.info("\n" + "="*60)
    logger.info("BUSCO Quality Summary:")
    logger.info("="*60)
    logger.info(f"\nCompleteness (mean): {report_df['completeness'].mean():.2f}%")
    logger.info(f"Completeness range: {report_df['completeness'].min():.2f}% - {report_df['completeness'].max():.2f}%")
    logger.info("="*60 + "\n")
    
    # Filter genomes
    passed, failed = filter_genomes(results, min_completeness=args.min_completeness)
    write_failed_genomes(failed)
    
    # Update metadata if requested
    if args.update_metadata and args.metadata:
        update_metadata(args.metadata, results)
    
    logger.info("\n✓ BUSCO quality assessment complete!")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Report file: busco_quality_report.tsv")
    if failed:
        logger.info(f"Failed genomes: genomes_failed_busco.tsv")


if __name__ == '__main__':
    main()
