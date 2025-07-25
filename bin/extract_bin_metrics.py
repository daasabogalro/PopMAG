#!/usr/bin/env python3

import argparse
import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union
import sys
import glob

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_contigs_to_bin(file_path: str) -> pd.DataFrame:
    """
    Read the contigs to bin file into a DataFrame.
    
    Args:
        file_path: Path to the contigs to bin file
        
    Returns:
        DataFrame containing contig to bin mapping
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        pd.errors.EmptyDataError: If the file is empty
    """
    try:
        return pd.read_csv(file_path, sep='\s+', header=None, names=['contig', 'bin'])
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except pd.errors.EmptyDataError:
        logger.error(f"File is empty: {file_path}")
        raise

def read_gene_metrics(file_path: str) -> pd.DataFrame:
    """
    Read the gene metrics file into a DataFrame.
    
    Args:
        file_path: Path to the gene metrics file
        
    Returns:
        DataFrame containing gene metrics
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        pd.errors.EmptyDataError: If the file is empty
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        # Extract sample name from file path and get only the part before "_to_"
        sample_name = Path(file_path).stem
        sample_id = sample_name.split('_to_')[0]
        df['sample_name'] = sample_id
        return df
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except pd.errors.EmptyDataError:
        logger.error(f"File is empty: {file_path}")
        raise

def get_bin_gene_metrics(contigs_to_bin_file: str, gene_metrics_files: List[str], bin_name: str) -> Optional[pd.DataFrame]:
    """
    Extract gene metrics for a specific bin from multiple gene metrics files.
    
    Args:
        contigs_to_bin_file: Path to the file containing contig to bin mapping
        gene_metrics_files: List of paths to gene metrics files
        bin_name: Name of the bin to extract metrics for
        
    Returns:
        DataFrame containing gene metrics for the specified bin or None if no contigs found
    """
    # Read input files
    contigs_df = read_contigs_to_bin(contigs_to_bin_file)
    
    # Filter contigs for the specified bin
    bin_contigs = contigs_df[contigs_df['bin'] == bin_name]['contig'].tolist()
    
    if not bin_contigs:
        logger.warning(f"No contigs found for bin {bin_name}")
        return None
    
    # Read and combine all gene metrics files
    all_metrics = []
    for metrics_file in gene_metrics_files:
        metrics_df = read_gene_metrics(metrics_file)
        all_metrics.append(metrics_df)
    
    combined_metrics = pd.concat(all_metrics, ignore_index=True)
    
    # Filter metrics based on bin contigs
    filtered_metrics = combined_metrics[combined_metrics['scaffold'].isin(bin_contigs)]
    return filtered_metrics

def calculate_sample_stats(df: pd.DataFrame) -> List[Dict[str, Union[str, int, float]]]:
    """
    Calculate statistics for each sample in the DataFrame.
    
    Args:
        df: DataFrame containing gene metrics
        
    Returns:
        List of dictionaries containing statistics for each sample
    """
    columns_to_sum = [
        'SNV_count', 'SNV_S_count', 'SNV_N_count', 
        'SNS_count', 'SNS_S_count', 'SNS_N_count'
    ]
    
    sample_stats = []
    for sample in df['sample_name'].unique():
        sample_df = df[df['sample_name'] == sample]
        stats = {col: sample_df[col].sum() for col in columns_to_sum}
        stats['sample_name'] = sample
        sample_stats.append(stats)
    
    return sample_stats

def process_multiple_bins(contigs_to_bin_file: str, gene_metrics_files: List[str], output_dir: str, sample_name: str) -> None:
    """
    Process multiple bins and create a combined summary file.
    
    Args:
        contigs_to_bin_file: Path to the file containing contig to bin mapping
        gene_metrics_files: List of paths to gene metrics files
        output_dir: Directory to save output files
    """
    # Read input files
    contigs_df = read_contigs_to_bin(contigs_to_bin_file)
    
    # Get unique bin names
    bin_names = contigs_df['bin'].unique()
    
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Initialize list to store summary statistics
    all_summaries = []
    
    # Process each bin
    for bin_name in bin_names:
        logger.info(f"Processing bin: {bin_name}")
        
        filtered_metrics = get_bin_gene_metrics(contigs_to_bin_file, gene_metrics_files, bin_name)
        
        if filtered_metrics is not None:
            # Calculate statistics for each sample
            sample_stats = calculate_sample_stats(filtered_metrics)
            
            # Add bin name to each sample's statistics
            for stats in sample_stats:
                stats['bin_name'] = bin_name
            
            all_summaries.extend(sample_stats)
            
            # Save individual bin metrics
            bin_output = output_path / f"{bin_name}_metrics.tsv"
            filtered_metrics = filtered_metrics.sort_values(['gene', 'sample_name'])
            filtered_metrics.to_csv(bin_output, sep='\t', index=False)
            logger.info(f"Saved {len(filtered_metrics)} gene metrics to {bin_output}")
    
    # Create combined summary file
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        # Reorder columns to put bin_name and sample_name first
        cols = ['bin_name', 'sample_name'] + [col for col in summary_df.columns if col not in ['bin_name', 'sample_name']]
        summary_df = summary_df[cols]
        
        combined_summary_file = output_path / f"SNVs_{sample_name}_summary.tsv"
        summary_df.to_csv(combined_summary_file, sep='\t', index=False)
        logger.info(f"Saved combined summary to {combined_summary_file}")
    else:
        logger.warning("No bins were processed successfully")

def main() -> None:
    """
    Main function to process gene metrics for bins.
    """
    parser = argparse.ArgumentParser(description="Extract bin metrics from inStrain gene_info files.")
    parser.add_argument("contigs_to_bins", help="Path to the contigs_to_bins.tsv file.")
    parser.add_argument("gene_metrics", nargs='+', help="Path to one or more gene_info.tsv files.")
    parser.add_argument("-d", "--output_dir", default=".", help="Output directory to save the metrics files.")
    parser.add_argument("-s", "--sample", required=True, help="Sample name for the summary file.")
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

    # Create the output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        process_multiple_bins(args.contigs_to_bins, args.gene_metrics, args.output_dir, args.sample)
    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=args.verbose)
        sys.exit(1)

if __name__ == "__main__":
    main()
