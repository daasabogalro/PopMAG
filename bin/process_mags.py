#!/usr/bin/env python3

from Bio import SeqIO
import argparse
from pathlib import Path
import sys

def process_mag_files(sample_id, mag_files):
    """
    Process multiple MAG FASTA files and write to a combined output file
    with MAG_name and Contig_name for each contig, and generate an individual
    definition file for each MAG with the same structure.

    Args:
        sample_id (str): Sample ID
        mag_files (list): List of paths to the input MAG FASTA files
    """
    combined_output_file = f"{sample_id}_combined_contigs.txt"
    
    with open(combined_output_file, 'w') as combined_out_f:
        for mag_file in mag_files:
            mag_name = Path(mag_file).stem  # Extract the base name of the MAG file
            individual_output_file = f"{sample_id}_{mag_name}_individual_contigs.txt"
            
            # Open the individual file for the current MAG with the same structure
            with open(individual_output_file, 'w') as mag_out_f:
                try:
                    # Iterate through the records (contigs) in the current MAG file
                    for record in SeqIO.parse(mag_file, "fasta"):
                        # Write the combined file with MAG and contig info
                        combined_line = f"{record.id}\t{mag_name}\n"
                        combined_out_f.write(combined_line)
                        
                        # Write the same line to the individual MAG file
                        mag_out_f.write(f"{record.id}\t{mag_name}\n")
                except Exception as e:
                    print(f"Error processing file {mag_file}: {str(e)}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Process multiple MAG FASTA files to extract MAG_name and Contig_name')
    parser.add_argument('sample_id', help='Sample ID')
    parser.add_argument('mag_files', nargs='+', help='MAG FASTA files')
    
    args = parser.parse_args()
    
    process_mag_files(args.sample_id, args.mag_files)

if __name__ == "__main__":
    main()