#!/usr/bin/env python

import pandas as pd
import requests
import numpy as np
from Bio import SeqIO

def get_markers(lineage, thr):
    # Get SNPs above threshold from Lapis and save to df
    url = f'https://lapis.cov-spectrum.org/open/v1/sample/nuc-mutations?pangoLineage={lineage}&downloadAsFile=false&dataFormat=json'
    print(f"Url is {url}")
    r = requests.get(url)
    result = r.json()
    df = pd.json_normalize(result['data'])

    try:                                                        # Removes SNPs below threshold
        df.drop(df[df['proportion'] < thr].index, inplace=True)
        df.drop(columns=['count'], inplace=True)
    except:                                                     # Unless no data returned from request
        print(f'No data available for {lineage}. Will be ignored.')
        return None

    # Get REF, POS and ALT from df 
    df[['REF', 'POS', 'ALT']] = df['mutation'].str.extract('(\D+)(\d+)(\D+)', expand=True)
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    df.sort_values(by=['POS'], ascending=True, inplace=True)

    return df

def consensus(df, reference):
    wuhan_seq = SeqIO.parse(open(reference),'fasta')            # Read Wuhan reference
    for fasta in wuhan_seq:                                     # Extract fasta sequence
        seq = list(str(fasta.seq))

    for _, row in df.iterrows():                                # Modify reference to create lineage consensus sequence
        pos = int(row['POS']) - 1
        alt = row['ALT']
        seq[pos] = alt

    return ''.join(seq)

def write_fasta(out_file, lineage, variant_fasta):
    with open(out_file + '.fasta', 'w') as of:
        ident = '>' + lineage + '\n'
        of.write(ident)                                         # Fasta identifier
        of.write(variant_fasta)                                 # Fasta sequence
    return None

def df_comparison(df_list, lineages):
    df1 = df_list[0][['mutation', 'proportion']]                # Save mutations of first lineage
    df2 = df_list[1][['mutation', 'proportion']]                # Save mutations of second lineage
    df = df1.merge(df2, on='mutation', how='outer', indicator=True)     # Merge all SNPs in one df

    df.replace({np.nan: '-',                                    # Rename values to ease understanding
                'left_only': lineages[0], 
                'right_only': lineages[1], 
                'both': 'Shared'}, inplace=True)

    df.rename(columns = {'proportion_x':'proportion_' + lineages[0], # Rename columns to ease understanding
                         'proportion_y':'proportion_' + lineages[1],
                         '_merge': 'Lineage'}, inplace = True)

    df[['REF', 'POS', 'ALT']] = df['mutation'].str.extract('(\D+)(\d+)(\D+)', expand=True)
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    df.sort_values(by=['POS'], ascending=True, inplace=True)
    col_list = ['mutation', 'REF', 'POS', 'ALT', 'proportion_BA.1', 'proportion_BA.2', 'Lineage']
    df = df[col_list]

    return df
