#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import logging
from pathlib import Path
from typing import Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def read_metrics_file(file_path: str) -> pd.DataFrame:
    """
    Read a metrics file into a DataFrame.
    
    Args:
        file_path: Path to the metrics file
        
    Returns:
        DataFrame containing metrics data
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        logger.info(f"Read metrics file: {file_path} with {len(df)} rows")
        return df
    except Exception as e:
        logger.error(f"Error reading metrics file {file_path}: {str(e)}")
        raise

def read_annotation_file(file_path: str) -> pd.DataFrame:
    """
    Read an annotation file into a DataFrame.
    
    Args:
        file_path: Path to the annotation file
        
    Returns:
        DataFrame containing annotation data
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        logger.info(f"Read annotation file: {file_path} with {len(df)} rows")
        return df
    except Exception as e:
        logger.error(f"Error reading annotation file {file_path}: {str(e)}")
        raise

def get_genome_name_from_metrics_file(file_path: str) -> str:
    """
    Extract genome name from metrics file path.
    
    Args:
        file_path: Path to the metrics file
        
    Returns:
        Genome name extracted from the file path
    """
    # Extract genome name from file path
    # Expected format: something like "genome_name_metrics.tsv"
    filename = os.path.basename(file_path)
    logger.debug(f"Extracting genome name from metrics file: {file_path}")
    logger.debug(f"Filename: {filename}")
    
    if filename.endswith("_metrics.tsv"):
        genome_name = filename.replace("_metrics.tsv", "")
        logger.debug(f"Found genome name from _metrics.tsv pattern: {genome_name}")
        return genome_name
    else:
        # Fallback: use filename without extension
        fallback_name = os.path.splitext(filename)[0]
        logger.debug(f"Using fallback genome name: {fallback_name}")
        return fallback_name

def merge_metrics_and_annotations(metrics_df: pd.DataFrame, annotations_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge metrics and annotations DataFrames.
    
    Args:
        metrics_df: DataFrame containing metrics data
        annotations_df: DataFrame containing annotation data
        
    Returns:
        Merged DataFrame
    """
    # Merge on gene ID
    # Try different possible column names for gene ID
    gene_cols_metrics = ['gene', 'gene_id', 'target']
    gene_cols_annotations = ['target', 'gene', 'gene_id']
    
    metrics_gene_col = None
    annotations_gene_col = None
    
    # Find the gene column in metrics
    for col in gene_cols_metrics:
        if col in metrics_df.columns:
            metrics_gene_col = col
            break
    
    # Find the gene column in annotations
    for col in gene_cols_annotations:
        if col in annotations_df.columns:
            annotations_gene_col = col
            break
    
    if metrics_gene_col is None:
        raise ValueError(f"Could not find gene column in metrics. Available columns: {list(metrics_df.columns)}")
    
    if annotations_gene_col is None:
        raise ValueError(f"Could not find gene column in annotations. Available columns: {list(annotations_df.columns)}")
    
    logger.info(f"Merging on metrics column '{metrics_gene_col}' and annotations column '{annotations_gene_col}'")
    
    # Merge the DataFrames
    merged = pd.merge(metrics_df, annotations_df, left_on=metrics_gene_col, right_on=annotations_gene_col, how='outer')
    
    # Remove specified columns if they exist
    columns_to_remove = ['scaffold', 'ORF_start', 'ORF_end', 'product_start', 'product_end']
    if annotations_gene_col != metrics_gene_col:
        columns_to_remove.append(annotations_gene_col)
    
    merged = merged.drop(columns=[col for col in columns_to_remove if col in merged.columns])
    
    return merged

def process_single_file_pair(bin_metrics_file: str, annotation_file: str, output_dir: str, genome_name: str = None) -> None:
    """
    Process a single pair of metrics and annotation files.
    
    Args:
        bin_metrics_file: Path to the metrics file
        annotation_file: Path to the annotation file
        output_dir: Directory to save output files
        genome_name: Optional genome name (if not provided, will be extracted from file)
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Processing single file pair:")
    logger.info(f"  Metrics file: {bin_metrics_file}")
    logger.info(f"  Annotation file: {annotation_file}")
    
    # Extract genome name if not provided
    if genome_name is None:
        genome_name = get_genome_name_from_metrics_file(bin_metrics_file)
    
    logger.info(f"  Genome name: {genome_name}")
    
    try:
        # Read the files
        metrics_df = read_metrics_file(bin_metrics_file)
        annotations_df = read_annotation_file(annotation_file)
        
        # Merge the data
        merged_df = merge_metrics_and_annotations(metrics_df, annotations_df)
        
        # Save the merged file
        output_file = output_path / f"{genome_name}_merged_metrics_with_annotation.tsv"
        merged_df.to_csv(output_file, sep='\t', index=False)
        
        logger.info(f"Successfully processed {genome_name}")
        logger.info(f"Merged file saved as {output_file}")
        logger.info(f"Columns in the final report: {list(merged_df.columns)}")
        logger.info(f"Total rows: {len(merged_df)}")
        
    except Exception as e:
        logger.error(f"Error processing {genome_name}: {str(e)}")
        raise

def main() -> None:
    """
    Main function to process and merge metrics with annotations.
    """
    parser = argparse.ArgumentParser(description='Merge inStrain metrics with MetaCerberus annotations')
    parser.add_argument('--bin_metrics', nargs='+', required=True, 
                       help='Path(s) to the bin metrics file(s)')
    parser.add_argument('--annotations', nargs='+', required=True,
                       help='Path(s) to the annotation file(s)')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for merged files')
    parser.add_argument('--genome_name', type=str,
                       help='Genome name for processing (optional)')
    
    args = parser.parse_args()
    
    try:
        # Validate input
        if len(args.bin_metrics) != 1 or len(args.annotations) != 1:
            raise ValueError("This script processes exactly one metrics file and one annotation file")
        
        process_single_file_pair(
            args.bin_metrics[0], 
            args.annotations[0], 
            args.output_dir,
            args.genome_name
        )
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        raise

if __name__ == "__main__":
    main() 