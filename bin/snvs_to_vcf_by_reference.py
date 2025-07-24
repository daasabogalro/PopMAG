import argparse
import os
import re
import sys
import pandas as pd
from collections import defaultdict
from pathlib import Path

def setup_logging(log_file=None):
    def log(msg):
        if log_file:
            with open(log_file, 'a') as lf:
                lf.write(f"{msg}\n")
        else:
            print(msg)
    return log

def parse_sample_info(filepath, sample_pattern=None):
    if sample_pattern:
        m = re.search(sample_pattern, os.path.basename(filepath))
    else:
        m = re.search(r"^(.+)_to_(.+)_reps\.IS_SNVs\.tsv$", os.path.basename(filepath))
    if m:
        return m.group(1), m.group(2)
    return None, None

def validate_input_files(snv_files, log):
    valid_files = []
    for f in snv_files:
        if not os.path.exists(f):
            log(f"Warning: File does not exist: {f}")
            continue
        if not os.access(f, os.R_OK):
            log(f"Warning: File not readable: {f}")
            continue
        valid_files.append(f)
    return valid_files

def process_snv_file(filepath, sample, log):
    try:
        df = pd.read_csv(filepath, sep='\t')
        df['sample'] = sample
        return df
    except Exception as e:
        log(f"Error reading {filepath}: {e}")
        return None

def write_vcf_header(vcf_file, samples, ref_name, contigs=None):
    vcf_file.write('##fileformat=VCFv4.2\n')
    vcf_file.write(f'##source=inStrain_SNV_converter\n')
    vcf_file.write(f'##reference={ref_name}\n')
    if contigs:
        for contig in contigs:
            vcf_file.write(f'##contig=<ID={contig}>\n')
    vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf_file.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">\n')
    vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')

def main():
    parser = argparse.ArgumentParser(
        description="Convert inStrain SNV files to VCFs grouped by reference genome set.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python snvs_to_vcf_by_reference.py --snv_files file1.tsv file2.tsv --out_dir ./vcf_output
        """
    )
    parser.add_argument('--snv_files', nargs='+', required=True,
                   help='List of SNV files to process (required)')
    parser.add_argument('--out_dir', type=str, default='./vcf_output', 
                       help='Directory to write VCF files (default: ./vcf_output)')
    parser.add_argument('--sample_pattern', type=str, default=None,
                       help='Custom regex pattern for parsing sample/reference from filename')
    parser.add_argument('--vcf_prefix', type=str, default='', 
                       help='Prefix for VCF filenames (default: "")')
    parser.add_argument('--vcf_suffix', type=str, default='', 
                       help='Suffix for VCF filenames (default: "")')
    parser.add_argument('--min_coverage', type=int, default=0, 
                       help='Minimum coverage (DP) to include a site (default: 0)')
    parser.add_argument('--min_var_freq', type=float, default=0, 
                       help='Minimum variant allele frequency to include a site (default: 0)')
    parser.add_argument('--include_nonvariant', action='store_true', 
                       help='Include non-variant sites (default: False)')
    parser.add_argument('--threads', type=int, default=1, 
                       help='Number of threads for parallel processing (default: 1)')
    parser.add_argument('--log_file', type=str, default=None, 
                       help='Log file path (default: None, print to stdout)')
    parser.add_argument('--verbose', action='store_true', 
                       help='Enable verbose logging')
    args = parser.parse_args()

    # Setup logging
    log = setup_logging(args.log_file)
    log(f"Starting SNV to VCF conversion")
    log(f"Input files: {args.snv_files}")
    snv_files = args.snv_files
    valid_files = validate_input_files(snv_files, log)
    if not valid_files:
        log(f"Error: No valid SNV files found")
        sys.exit(1)

    log(f"Output directory: {args.out_dir}")
    try:
        os.makedirs(args.out_dir, exist_ok=True)
    except Exception as e:
        log(f"Error creating output directory: {e}")
        sys.exit(1)

    log(f"Found {len(valid_files)} valid SNV files")

    # Group SNV files by reference genome set
    ref_to_files = defaultdict(list)
    for f in valid_files:
        sample, ref = parse_sample_info(f, args.sample_pattern)
        if sample and ref:
            ref_to_files[ref].append((sample, f))
        else:
            log(f"Warning: Could not parse sample/reference from {f}")

    if not ref_to_files:
        log(f"Error: No valid sample/reference pairs found")
        sys.exit(1)
    
    # Process each reference set
    processed_refs = 0
    for ref, sample_files in ref_to_files.items():
        log(f"\nProcessing reference set: {ref} ({len(sample_files)} samples)")
        
        # Read and process SNV files
        dfs = {}
        for sample, f in sample_files:
            if args.verbose:
                log(f"  Reading {sample} from {f}")
            df = process_snv_file(f, sample, log)
            if df is not None:
                dfs[sample] = df
        
        if not dfs:
            log(f"  No valid data for reference {ref}, skipping")
            continue
        
        # Merge SNV tables
        merged = None
        for sample, df in dfs.items():
            allele_cols = ['A', 'C', 'T', 'G']
            rename_dict = {col: f'{col}_{sample}' for col in allele_cols}
            df = df.rename(columns=rename_dict)
            
            # Rename ref_base and var_base to avoid conflicts
            if 'ref_base' in df.columns:
                df = df.rename(columns={'ref_base': f'ref_base_{sample}'})
            if 'var_base' in df.columns:
                df = df.rename(columns={'var_base': f'var_base_{sample}'})
            
            # Select columns for merging
            merge_cols = ['scaffold', 'position']
            allele_cols_renamed = [f'{b}_{sample}' for b in allele_cols]
            ref_var_cols = []
            if f'ref_base_{sample}' in df.columns:
                ref_var_cols.append(f'ref_base_{sample}')
            if f'var_base_{sample}' in df.columns:
                ref_var_cols.append(f'var_base_{sample}')
            
            df_to_merge = df[merge_cols + allele_cols_renamed + ref_var_cols]
            
            if merged is None:
                merged = df_to_merge.copy()
            else:
                merged = pd.merge(merged, df_to_merge, on=merge_cols, how='outer')
        
        # Fill missing allele counts with 0
        for sample in dfs.keys():
            for base in ['A', 'C', 'T', 'G']:
                col = f'{base}_{sample}'
                if col not in merged.columns:
                    merged[col] = 0
                merged[col] = merged[col].fillna(0).astype(int)
        
        # Handle ref_base and var_base consolidation
        ref_base_cols = [col for col in merged.columns if col.startswith('ref_base_')]
        var_base_cols = [col for col in merged.columns if col.startswith('var_base_')]
        
        merged['ref_base'] = merged[ref_base_cols].bfill(axis=1).iloc[:, 0] if ref_base_cols else '.'
        merged['var_base'] = merged[var_base_cols].bfill(axis=1).iloc[:, 0] if var_base_cols else '.'
        
        # Convert 0-based to 1-based positions for VCF
        merged['POS'] = merged['position'] + 1
        samples = list(dfs.keys())
        
        # Sort merged DataFrame by scaffold and POS
        merged_sorted = merged.sort_values(['scaffold', 'POS'])
        
        # Write VCF file
        vcf_name = f"{args.vcf_prefix}{ref}{args.vcf_suffix}.vcf"
        vcf_path = os.path.join(args.out_dir, vcf_name)
        log(f"  Writing VCF: {vcf_path} with {len(samples)} samples")
        
        # Collect unique contig names
        contigs = list(merged['scaffold'].dropna().unique())
        
        try:
            with open(vcf_path, 'w') as vcf:
                write_vcf_header(vcf, samples, ref, contigs=contigs)
                
                sites_written = 0
                for idx, row in merged_sorted.iterrows():
                    chrom = row['scaffold']
                    pos = row['POS']
                    
                    # Get consolidated ref_base and var_base
                    ref_base = row.get('ref_base', '.')
                    alt_base = row.get('var_base', '.')
                    
                    # Skip sites with undetermined alleles
                    if pd.isna(ref_base) or pd.isna(alt_base) or ref_base == '.' or alt_base == '.':
                        continue
                    
                    # Handle SNSs - find most common non-ref allele if var_base equals ref_base
                    if alt_base == ref_base:
                        allele_totals = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
                        for sample in samples:
                            for base in ['A', 'C', 'T', 'G']:
                                if base != ref_base:
                                    allele_totals[base] += row.get(f'{base}_{sample}', 0)
                        
                        max_count = 0
                        most_common_alt = None
                        for base, count in allele_totals.items():
                            if count > max_count:
                                max_count = count
                                most_common_alt = base
                        
                        if most_common_alt and max_count > 0:
                            alt_base = most_common_alt
                        else:
                            continue
                    
                    # Check for variation
                    has_variation = False
                    for sample in samples:
                        total_reads = sum(row.get(f'{b}_{sample}', 0) for b in ['A', 'C', 'T', 'G'])
                        if total_reads > 0:
                            allele_counts = [row.get(f'{b}_{sample}', 0) for b in ['A', 'C', 'T', 'G']]
                            non_zero_alleles = [count for count in allele_counts if count > 0]
                            if len(non_zero_alleles) > 1:
                                has_variation = True
                                break
                    
                    # Apply filters
                    if not args.include_nonvariant and not has_variation:
                        continue
                    
                    # Coverage filter
                    if args.min_coverage > 0:
                        total_coverage = 0
                        for sample in samples:
                            for base in ['A', 'C', 'T', 'G']:
                                total_coverage += row.get(f'{base}_{sample}', 0)
                        if total_coverage < args.min_coverage:
                            continue
                    
                    # Variant frequency filter
                    if args.min_var_freq > 0:
                        max_var_freq = 0
                        for sample in samples:
                            total_reads = sum(row.get(f'{b}_{sample}', 0) for b in ['A', 'C', 'T', 'G'])
                            if total_reads > 0:
                                ref_count = row.get(f'{ref_base}_{sample}', 0)
                                var_count = total_reads - ref_count
                                if var_count > 0:
                                    var_freq = var_count / total_reads
                                    max_var_freq = max(max_var_freq, var_freq)
                        if max_var_freq < args.min_var_freq:
                            continue
                    
                    # Write VCF line
                    format_field = 'GT:AD'
                    sample_fields = []
                    for sample in samples:
                        ref_count = row.get(f'{ref_base}_{sample}', 0)
                        alt_count = row.get(f'{alt_base}_{sample}', 0)
                        ad = f"{ref_count},{alt_count}"
                        gt = './.'  # metagenomic: unknown genotype
                        sample_fields.append(f"{gt}:{ad}")
                    
                    vcf.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\t{format_field}\t" + '\t'.join(sample_fields) + '\n')
                    sites_written += 1
                
                log(f"  Wrote {sites_written} sites to {vcf_path}")
                processed_refs += 1
                
        except Exception as e:
            log(f"Error writing VCF file {vcf_path}: {e}")
            continue
    
    if processed_refs == 0:
        log("Error: No VCF files were successfully created")
        sys.exit(1)
    
    log(f"\nSuccessfully processed {processed_refs} reference sets")
    log("All VCFs written.")

if __name__ == "__main__":
    main() 