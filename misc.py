#!/usr/bin/env python

import pandas as pd
import requests
from Bio import SeqIO

def get_markers(lineage, thr):
    # Get SNPs above threshold from Lapis and save to df
    url = f'https://lapis.cov-spectrum.org/open/v1/sample/nuc-mutations?pangoLineage={lineage}&downloadAsFile=false&dataFormat=json'
    print(f"Url is {url}")
    r = requests.get(url)
    result = r.json()
    df = pd.json_normalize(result['data'])
    # Removes SNPs below threshold
    try:
        df.drop(df[df['proportion'] < thr].index, inplace=True)
        df.drop(columns=['count'], inplace=True)
    # Unless no data returned from request
    except:
        print(f'No data available for {lineage}. Will be ignored.')
        return None

    # Get REF, POS and ALT from df 
    df[['REF', 'POS', 'ALT']] = df['mutation'].str.extract('(\D+)(\d+)(\D+)', expand=True)
    df.sort_values(by=['POS'], ascending=False, inplace=True)

    return df

def consensus(df, reference):
    # Read Wuhan reference
    wuhan_seq = SeqIO.parse(open(reference),'fasta')
    # Extract fasta sequence
    for fasta in wuhan_seq:
        seq = list(str(fasta.seq))
    # Modify reference to create lineage sequence
    for _, row in df.iterrows():
        pos = int(row['POS']) - 1
        alt = row['ALT']
        seq[pos] = alt

    return ''.join(seq)

def write_fasta(out_file, lineage, variant_fasta):
    with open(out_file + '.fasta', 'w') as of:
        ident = '>' + lineage + '\n'
        of.write(ident)
        of.write(variant_fasta)
    return None