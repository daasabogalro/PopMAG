#!/usr/bin/env python3

import argparse
import pandas as pd
import logging
from pathlib import Path
from typing import List, Dict
import sys
import glob

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_genome_info(file_path: str) -> pd.DataFrame:
    """
    Read a genome info file and add sample information.
    
    Args:
        file_path: Path to the genome info file
        
    Returns:
        DataFrame containing genome information with sample name
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        # Extract sample name from file path (part before _to_)
        sample_name = Path(file_path).stem.split('_to_')[0]
        df['sample'] = sample_name
        return df
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {str(e)}")
        raise

def get_reference_sample(file_path: str) -> str:
    """
    Extract the reference sample name from the file path.
    
    Args:
        file_path: Path to the genome info file
        
    Returns:
        Reference sample name (part between _to_ and _reps)
    """
    try:
        # Extract the part between _to_ and _reps
        filename = Path(file_path).stem
        reference = filename.split('_to_')[1].split('_reps')[0]
        return reference
    except Exception as e:
        logger.error(f"Error extracting reference sample from {file_path}: {str(e)}")
        raise

def combine_genome_info(files: List[str], output_file: str) -> None:
    """
    Combine genome information from multiple files into a single report.
    
    Args:
        files: List of genome info file paths
        output_file: Path to save the combined report
    """
    try:
        # Get reference sample from first file
        reference_sample = get_reference_sample(files[0])
        
        # Modify output filename to include reference sample
        output_path = Path(output_file)
        output_file = str(output_path.parent / f"{output_path.stem}_{reference_sample}{output_path.suffix}")
        
        # Read and combine all files
        dfs = []
        for file in files:
            logger.info(f"Processing file: {file}")
            df = read_genome_info(file)
            dfs.append(df)
        
        # Combine all DataFrames
        combined_df = pd.concat(dfs, ignore_index=True)
        
        # Reorder columns to put sample and genome first
        cols = ['sample', 'genome'] + [col for col in combined_df.columns if col not in ['sample', 'genome']]
        combined_df = combined_df[cols]
        
        # Sort by genome first, then by sample
        combined_df = combined_df.sort_values(['genome', 'sample'])
        
        # Save combined report
        combined_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved combined report to {output_file}")
        
        # Print summary statistics
        logger.info("\nSummary of samples and genomes:")
        sample_counts = combined_df.groupby('sample')['genome'].count()
        logger.info("\nNumber of genomes per sample:")
        for sample, count in sample_counts.items():
            logger.info(f"{sample}: {count} genomes")
        
        # Calculate and print average metrics per sample
        metrics = ['coverage', 'breadth', 'nucl_diversity', 'length', 'SNV_count', 'SNS_count']
        logger.info("\nAverage metrics per sample:")
        for metric in metrics:
            if metric in combined_df.columns:
                avg_by_sample = combined_df.groupby('sample')[metric].mean()
                logger.info(f"\nAverage {metric} per sample:")
                for sample, value in avg_by_sample.items():
                    logger.info(f"{sample}: {value:.4f}")
        
    except Exception as e:
        logger.error(f"Error combining genome info: {str(e)}")
        sys.exit(1)

def main() -> None:
    """
    Main function to combine genome information from multiple files.
    """
    parser = argparse.ArgumentParser(description='Combine genome information from multiple samples')
    parser.add_argument('files', nargs='+', help='Path(s) to genome info file(s)')
    parser.add_argument('-o', '--output', required=True, help='Output file path for combined report')
    
    args = parser.parse_args()
    
    combine_genome_info(args.files, args.output)

if __name__ == "__main__":
    main() 