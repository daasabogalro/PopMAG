#!/usr/bin/env python3

import csv
import sys

def transform_report(input_file, output_file, extension):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile)
        
        writer.writerow(['genome', 'completeness', 'contamination'])
        
        next(reader)  # Skip header
        
        for row in reader:
            genome = f"{row[0]}.{extension}"
            completeness = float(row[1])
            contamination = float(row[2])
            
            writer.writerow([genome, f"{completeness:.2f}", f"{contamination:.2f}"])

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extension = sys.argv[3]
    
    transform_report(input_file, output_file, extension)
