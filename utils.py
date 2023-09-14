#!/usr/bin/env python

import pandas as pd
import requests
import numpy as np
from Bio import SeqIO
import warnings
from pandas.errors import SettingWithCopyWarning
import simple_colors

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)

def get_markers(lineage, thr):
    # Get SNPs above threshold from Lapis and save to df
    url = f'https://lapis.cov-spectrum.org/open/v1/sample/nuc-mutations?pangoLineage={lineage}&downloadAsFile=false&dataFormat=json'
    print(f"Getting {lineage} SNPs", end=' ')
    r = requests.get(url)
    result = r.json()
    df = pd.json_normalize(result['data'])

    try:                                                        # Removes SNPs below threshold
        df.drop(df[df['proportion'] < thr].index, inplace=True)
        df.drop(columns=['count'], inplace=True)
        print(simple_colors.green(f'--> OK','bold'))
    except:                                                     # Unless no data returned from request
        print(simple_colors.red(f'--> No data available.','bold'))
        return None

    # Get REF, POS and ALT from df 
    df[['REF', 'POS', 'ALT']] = df['mutation'].str.extract('(\D+)(\d+)(\D+)', expand=True)
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    df.sort_values(by=['POS'], ascending=True, inplace=True)
    df.name = lineage

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

def df_comparison(df_list, lineages, thr):

    if len(df_list) < 2:
        print(simple_colors.red(f"At least 2 DataFrames are required for comparison. Comparison could not be performed.\nPlease, check output.", 'bold'))
        exit()

    in_list = [f'in_{lineages[0]}']
    proportion_list = [f'proportion_{lineages[0]}']

    merged_df = df_list[0][['mutation', 'proportion']]
    conditions = [[merged_df['proportion'] < thr, merged_df['proportion'] >= thr], [0, 1]]
    
    merged_df = merged_df.copy()
    merged_df = df_list[0][['mutation', 'proportion']]
    merged_df[f'in_{lineages[0]}'] = np.select(conditions[0], conditions[1])
    merged_df.rename(columns={'proportion': f'proportion_{lineages[0]}'}, inplace=True)

    # Merge the remaining DataFrames
    for i in range(1, len(df_list)):
        df_i = df_list[i][['mutation', 'proportion']]
        conditions = [[df_i['proportion'] < thr, df_i['proportion'] >= thr], [0, 1]]
        df_i[f'in_{lineages[i]}'] = np.select(conditions[0], conditions[1])
        in_list.append(f'in_{lineages[i]}')
        df_i.rename(columns={'proportion': f'proportion_{lineages[i]}'}, inplace=True)
        proportion_list.append(f'proportion_{lineages[i]}')
        # Merge on 'mutation' column using outer join
        merged_df = pd.merge(merged_df, df_i, on='mutation', how='outer', indicator=False)

    # Extract REF, POS, and ALT columns
    merged_df[['REF', 'POS', 'ALT']] = merged_df['mutation'].str.extract('(\D+)(\d+)(\D+)', expand=True)
    merged_df['POS'] = pd.to_numeric(merged_df['POS'], errors='coerce')

    # Sort by 'POS' column
    merged_df.sort_values(by=['POS'], ascending=True, inplace=True)
    merged_df.replace({np.nan: 0}, inplace=True)

     # Sort columns
    cols = ['mutation', 'REF', 'POS', 'ALT'] + proportion_list + in_list
    merged_df = merged_df[cols]
    merged_df = merged_df.loc[~(merged_df[in_list].eq(0).all(axis=1))]

    return merged_df
