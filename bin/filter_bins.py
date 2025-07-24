#!/usr/bin/env python3

import csv
import os
import shutil
import sys
import re

def filter_bins(checkm2_output_files, transformed_report, meta_id, min_completeness, max_contamination):
    os.mkdir(f"{meta_id}_filtered_bins")

    with open(transformed_report, 'r') as f:
        reader = csv.DictReader(f)
        with open(f"filtered_report_{meta_id}.csv", 'w', newline='') as out_f:
            fieldnames = ['genome', 'completeness', 'contamination']
            writer = csv.DictWriter(out_f, fieldnames=fieldnames)
            writer.writeheader()
            
            for row in reader:
                completeness = float(row['completeness'])
                contamination = float(row['contamination'])
                
                if completeness >= min_completeness and contamination <= max_contamination:
                    writer.writerow(row)
                    genome_name = row['genome']
                    matching_files = [f for f in checkm2_output_files if genome_name in f]
                    if matching_files:
                        genome_file = matching_files[0]
                        dest_path = os.path.join(f"{meta_id}_filtered_bins", os.path.basename(genome_file))
                        shutil.copy(genome_file, dest_path)
                        prefix = os.path.splitext(os.path.basename(genome_file))[0]
                        # Sanitize prefix: replace non-alphanumeric (except _) with _
                        prefix = re.sub(r'[^A-Za-z0-9_]', '_', prefix)
                        prefix_contig_headers(dest_path, prefix)

def prefix_contig_headers(fasta_path, prefix):
    temp_path = fasta_path + ".tmp"
    contig_num = 1
    with open(fasta_path, 'r') as infile, open(temp_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                new_header = f">{prefix}_contig{contig_num:09d}\n"
                outfile.write(new_header)
                contig_num += 1
            else:
                outfile.write(line)
    os.replace(temp_path, fasta_path)

if __name__ == "__main__":
    checkm2_output_files = sys.argv[1].split()
    transformed_report = sys.argv[2]
    meta_id = sys.argv[3]
    min_completeness = float(sys.argv[4])
    max_contamination = float(sys.argv[5])
    
    filter_bins(checkm2_output_files, transformed_report, meta_id, min_completeness, max_contamination)