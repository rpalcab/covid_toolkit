#!/usr/bin/env python

# Third party imports
import argparse
import os
import warnings
import simple_colors

# Local imports
from utils import get_markers, df_comparison

def get_arguments():

        parser = argparse.ArgumentParser(
            prog='lineage_comparison.py', description='Script to ease COVID lineages comparison')

        input_group = parser.add_argument_group('Input', 'Input parameters')
        input_group.add_argument('-l', '--lin', metavar="Lineages", action='store', required=True, help='REQUIRED. Determine lineages using quotes. Two lineages are required')
        input_group.add_argument('-t', '--thr', metavar="Threshold", action='store', type=float, default=0.9, help='Determine SNP frequency threshold. Min=0, max=1. Default=0.9')

        output_group = parser.add_argument_group(
            'Output', 'Required parameter to output results')
        output_group.add_argument('-o', '--out', metavar="Out directory", action='store', type=str, required=True, default='.', help='REQUIRED. Output directory to extract all results')

        arguments = parser.parse_args()

        return arguments

def main():
        # Get arguments from command line
        args = get_arguments()
        lineages = [i for i in args.lin.split(' ')]
        if len(lineages) >= 2:
            pass
        else:
            exit('Please, introduce two or more COVID lineages')
        thr = args.thr
        out_path = os.path.abspath(args.out)

        df_list = []

        # Get marker SNPs of two lineages
        for i, lineage in enumerate(lineages):
            df_list.append(get_markers(lineage, 0))
            if df_list[i] is None:                                  # Warn if no data of a lineage
                warnings.warn(f'Lineage {lineage} is not valid.')

        # Lineages comparison
        df_list = [df for df in df_list if df is not None]
        avail_lins = [df.name for df in df_list]
        comp_df = df_comparison(df_list, avail_lins, thr)
        out_file = os.path.join(out_path, '_'.join([i for i in avail_lins]))
        comp_df.to_csv(out_file + '.csv', index=False)
        
        print(simple_colors.green('Comparison finished!', 'bold'))

if __name__ == "__main__":
        main()
