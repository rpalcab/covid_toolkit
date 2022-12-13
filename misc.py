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
    df.drop(df[df['proportion'] < thr].index, inplace=True)
    df.drop(columns=['proportion', 'count'], inplace=True)

    # Get REF, POS and ALT from df 
    df[['REF', 'POS', 'ALT']] = df['mutation'].str.extract('(\D+)(\d+)(\D+)', expand=True)

    return df

def consensus(df):
    # Read Wuhan reference
    wuhan_seq = SeqIO.parse(open('/home/laura/iisgm/NC_045512.2.fasta'),'fasta')
    # Extract fasta sequence
    for fasta in wuhan_seq:
        seq = list(str(fasta.seq))
    # Modify reference to create lineage sequence
    for _, row in df.iterrows():
        pos = int(row['POS']) - 1
        alt = row['ALT']
        seq[pos] = alt

    return ''.join(seq)