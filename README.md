# COVID TOOLKIT
A set of Python scripts to extract relevant genomic information from COVID lineages. Uses Lapis to collect data of SNPs present in GISAID sequence database.

## Programs
### extract_snps.py
Recovers marker SNPs of one or more lineages. The output is a .csv file indicating the SNPs and their frequency in that lineage based on GISAID samples. Also, it generates a fasta file with the consensus sequence of the lineages.

The user can modify the minimum threshold frequency of SNPs to be considered (0.9 as default value). One or several lineages can be used as input, as well as a group of sublineages (e. g.: BA.2*). The COVID reference used to create the lineages consensus can also be changed (NC_045512.2 by default).

**Usage:**

```python extract_snps.py -l "BA.2.1* BQ.1" [-t 0.75] [-r /path/to/reference/NC_045512.2.fasta] -o .```

### lineage_comparison.py
Gets differential SNPs between two SARS-CoV19 lineages or sublineages. It creates a .csv file of the marker SNPs of both lineages, indicating to which one of them (or both) they belong.

Similarly to extract_snps.py, it can accept specific lineages or a set of sublineages (marked with an asterisk), as well as change the minimum threshold frequency of the SNPs.

**Usage:**

```python lineage_comparison.py -l "BA.2.1* BQ.1" -t 0.75 -o .```

## Next steps
- [ ] Merge with lineage_tracking
 
