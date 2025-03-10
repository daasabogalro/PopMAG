#!/usr/bin/env python3

from Bio import SeqIO
import argparse
from pathlib import Path
import sys

def process_mag_files(sample_id, mag_files):
    """
    Process multiple MAG FASTA files and write to a combined output file
    with MAG_name and Contig_name for each contig.

    Args:
        sample_id (str): Sample ID
        mag_files (list): List of paths to the input MAG FASTA files
    """
    output_file = f"{sample_id}_combined_contigs.txt"
    
    with open(output_file, 'w') as out_f:
        for mag_file in mag_files:
            mag_name = Path(mag_file).stem
            try:
                for record in SeqIO.parse(mag_file, "fasta"):
                    line = f"{record.id}\t{mag_name}\n"
                    out_f.write(line)
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