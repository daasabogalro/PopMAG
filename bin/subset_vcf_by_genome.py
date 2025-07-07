#!/usr/bin/env python3
"""
Subset a single VCF file by genome using a corresponding contigs2bin.tsv file.
For use in Nextflow workflows: processes one VCF and one contigs2bin.tsv at a time.
"""

import argparse
import pandas as pd
import os
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Subset a VCF file by genome using a contigs2bin.tsv file")
    parser.add_argument('--vcf_file', type=str, required=True, help='Input VCF file')
    parser.add_argument('--contigs2bin_file', type=str, required=True, help='contigs2bin.tsv file')
    parser.add_argument('--out_dir', type=str, required=True, help='Output directory for genome-specific VCFs')
    parser.add_argument('--prefix', type=str, default='', help='Prefix for output VCF filenames')
    parser.add_argument('--suffix', type=str, default='', help='Suffix for output VCF filenames')
    parser.add_argument('--include_unassigned', action='store_true', help='Include contigs not in contigs2bin.tsv in an "unassigned" genome VCF')
    return parser.parse_args()

def read_contigs2bin(contigs2bin_file):
    """Read contigs2bin.tsv and return mapping of contig to genome"""
    try:
        df = pd.read_csv(contigs2bin_file, sep='\t', header=None, names=['contig', 'genome'])
    except:
        try:
            df = pd.read_csv(contigs2bin_file, sep=',', header=None, names=['contig', 'genome'])
        except:
            df = pd.read_csv(contigs2bin_file, sep='\t')
            if len(df.columns) >= 2:
                df = df.iloc[:, [0, 1]]
                df.columns = ['contig', 'genome']
            else:
                raise ValueError("contigs2bin.tsv must have at least 2 columns")
    return dict(zip(df['contig'], df['genome']))

def convert_ad_to_ro_ao(vcf_line):
    parts = vcf_line.strip().split('\t')
    if len(parts) < 9:
        return vcf_line
    format_field = parts[8]
    sample_fields = parts[9:]
    new_format = format_field.replace('AD', 'RO:AO')
    new_sample_fields = []
    for sample_field in sample_fields:
        if sample_field == '.':
            new_sample_fields.append('.')
            continue
        sample_parts = sample_field.split(':')
        format_parts = format_field.split(':')
        try:
            ad_index = format_parts.index('AD')
            ad_value = sample_parts[ad_index]
        except (ValueError, IndexError):
            new_sample_fields.append(sample_field)
            continue
        if ad_value == '.':
            new_sample_fields.append(sample_field.replace(ad_value, '.:.'))
            continue
        try:
            ad_counts = ad_value.split(',')
            if len(ad_counts) >= 2:
                ref_count = ad_counts[0]
                alt_count = ad_counts[1]
                new_sample_parts = list(sample_parts)
                new_sample_parts[ad_index] = f"{ref_count}:{alt_count}"
                new_sample_fields.append(':'.join(new_sample_parts))
            else:
                new_sample_fields.append(sample_field)
        except:
            new_sample_fields.append(sample_field)
    new_parts = parts[:8] + [new_format] + new_sample_fields
    return '\t'.join(new_parts) + '\n'

def subset_vcf_by_genome(vcf_file, contig_to_genome, out_dir, prefix='', suffix='', include_unassigned=False):
    os.makedirs(out_dir, exist_ok=True)
    genome_variants = defaultdict(list)
    skipped_contigs = set()
    total_variants = 0
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            total_variants += 1
            chrom = parts[0]
            genome = contig_to_genome.get(chrom)
            if genome is None:
                skipped_contigs.add(chrom)
                if include_unassigned:
                    genome = "unassigned"
                else:
                    continue
            converted_line = convert_ad_to_ro_ao(line)
            genome_variants[genome].append(converted_line)
    header_lines = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##FORMAT=<ID=AD'):
                header_lines.append('##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">\n')
                header_lines.append('##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">\n')
            elif line.startswith('#'):
                header_lines.append(line)
            else:
                break
    genome_counts = defaultdict(int)
    for genome, variants in genome_variants.items():
        if not variants:
            continue
        vcf_filename = f"{prefix}{genome}{suffix}.vcf"
        vcf_path = os.path.join(out_dir, vcf_filename)
        with open(vcf_path, 'w') as current_vcf:
            for header_line in header_lines:
                current_vcf.write(header_line)
            for variant_line in variants:
                current_vcf.write(variant_line)
                genome_counts[genome] += 1

def main():
    args = parse_arguments()
    contig_to_genome = read_contigs2bin(args.contigs2bin_file)
    subset_vcf_by_genome(
        args.vcf_file,
        contig_to_genome,
        args.out_dir,
        args.prefix,
        args.suffix,
        args.include_unassigned
    )

if __name__ == "__main__":
    main() 