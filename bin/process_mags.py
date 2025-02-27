#!/usr/bin/env python3

from Bio import SeqIO
import argparse
from pathlib import Path

def process_mag_file(input_file, output_file):
    """
    Process a single MAG FASTA file and append to the output file
    with MAG_name and Contig_name for each contig.

    Args:
        input_file (str): Path to the input MAG FASTA file
        output_file (str): Path to output file
    """
    mag_name = Path(input_file).stem
    
    with open(output_file, 'a') as out_f:
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                # Write mag_name and contig name to output file
                out_f.write(f"{record.id}\t{mag_name}\n")
        except Exception as e:
            print(f"Error processing file {input_file}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Process MAG FASTA file to extract MAG_name and Contig_name')
    parser.add_argument('input_file', help='MAG FASTA file')
    parser.add_argument('output_file', help='Output file path')
    
    args = parser.parse_args()
    
#    # If the output file doesn't exist, create it with a header
#    if not os.path.exists(args.output_file):
#        with open(args.output_file, 'w') as out_f:
#            out_f.write()
    
    process_mag_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()