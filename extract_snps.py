#!/usr/bin/env python

# Third party imports
import argparse
import os
import simple_colors

# Local imports
from utils import get_markers, consensus, write_fasta

def get_arguments():

        parser = argparse.ArgumentParser(
            prog='extract_snps.py', description='Toolkit to ease COVID analysis and comparison of lineages and sublineages')

        input_group = parser.add_argument_group('Input', 'Input parameters')
        input_group.add_argument('-l', '--lin', metavar="Lineages", action='store', required=True, help='REQUIRED. Determine lineage(s). If two or more, use quotes')
        input_group.add_argument('-t', '--thr', metavar="Threshold", action='store', type=float, default=0.9, help='Determine threshold. Min=0, max=1, def=0.9')
        input_group.add_argument('-r', '--ref', metavar="Reference", action='store', type=str, default=os.path.join(os.path.realpath(os.path.dirname(__file__)), 'NC_045512.2.fasta'), help='COVID reference. Default: NC_045512.2')

        output_group = parser.add_argument_group(
            'Output', 'Required parameter to output results')
        output_group.add_argument('-o', '--out', metavar="Out directory", action='store', type=str, required=True, default='.', help='REQUIRED. Output directory to extract all results')

        arguments = parser.parse_args()

        return arguments

# Get arguments from comment line
args = get_arguments()
lineages = [i for i in args.lin.split(' ')]
thr = args.thr
out_path = os.path.abspath(args.out)
ref = args.ref

# Check if reference file exists
if os.path.exists(ref):
    pass
else:
    print(simple_colors.red(f'Error. Reference file not found. Please check {ref}', 'bold'))
    exit()

# Extract SNPs of every lineage
for lineage in lineages:
    df = get_markers(lineage, thr)
    if df is None:                                  # Next iteration if no SNP data
        continue
    variant_fasta = consensus(df, ref)

    out_file = os.path.join(out_path, lineage)
    df.to_csv(out_file + '.csv', index=False)       # Outuput SNP table to new_file
    write_fasta(out_file, lineage, variant_fasta)   # Output sequence to new file

print(simple_colors.green('Execution finished!', 'bold'))
