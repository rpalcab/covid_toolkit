#!/usr/bin/env python

# Third party imports
import argparse
import os

# Local imports
from misc import get_markers, consensus

def get_arguments():

        parser = argparse.ArgumentParser(
            prog='extract_snps.py', description='Toolkit to ease COVID analysis and comparison of lineages and sublineages')

        input_group = parser.add_argument_group('Input', 'Input parameters')
        input_group.add_argument('-l', '--lin', metavar="Lineages", action='store', required=True, help='REQUIRED. Determine lineage(s). If two or more, use quotes')
        input_group.add_argument('-t', '--thr', metavar="Threshold", action='store', type=float, default=0.9, help='Determine threshold. Min=0, max=1')

        output_group = parser.add_argument_group(
            'Output', 'Required parameter to output results')
        output_group.add_argument('-o', '--out', type=str, required=True, default='.', help='REQUIRED. Output directory to extract all results')

        arguments = parser.parse_args()

        return arguments

args = get_arguments()

lineages = [i for i in args.lin.split(' ')]
thr = args.thr
out_path = os.path.abspath(args.out)

for lineage in lineages:
    df = get_markers(lineage, thr)
    variant_fasta = consensus(df)

    out_file = os.path.join(out_path, lineage)
    # Outuput SNP table to new_file
    df.to_csv(out_file + '.csv', index=False)
    # Output sequence to new file 
    
    with open(out_file + '.fasta', 'w') as of:
        ident = '>' + lineage + '\n'
        of.write(ident)
        of.write(variant_fasta)
