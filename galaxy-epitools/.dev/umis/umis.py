#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import os
import regex as re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


# ················································································· #

def umis(seq: str,
         umi_5_pattern: str,
         umi_3_pattern: str,
         n_errors: int=0) -> tuple:
    """
    Takes a sequence, tries to find UMIs within a certain margin of error,
    and returns those UMIs.

    :param seq: A sequence in which to search for the UMI in.
    :param umi_5_pattern: The pattern for a 5'-end UMI.
    :param umi_3_pattern: The pattern for a 3'-end UMI
    :param n_errors: How many errors to accept when searching for either UMI.
    :param not_found: What to return, when a UMI is not found.
    :return: A tuple of (5' UMI, 3' UMI). If UMIs are not found, value will be ``None``.
    """
    umi_5_pattern = '(' + umi_5_pattern + '){s<=' + str(n_errors) + '}'
    umi_3_pattern = '(' + umi_3_pattern + '){s<=' + str(n_errors) + '}'
    umi_matches = (re.search(pat, seq) for pat in [umi_5_pattern, umi_3_pattern])
    return tuple(match.group() if match is not None else None for match in umi_matches)



def extract_umis(input_fastq: str=None,
                 input_fasta: str=None,
                 output_tsv: str=None,
                 umi_5_pattern: str=r'T{3}[ACG]{4}T{2}[ACG]{4}T{2}[ACG]{4}T{2}[ACG]{4}T{3}',
                 umi_3_pattern: str=r'A{3}[TCG]{4}A{2}[TCG]{4}A{2}[TCG]{4}A{2}[TCG]{4}A{3}',
                 n_errors: int=2,
                 not_found: str='N/A',
                 unique_by_revcomp: bool=True):
    """
    Extracts 5' and 3' UMIs from sequencing reads in a FASTQ or FASTA file and saves them in a TSV file.

    :param input_fastq: Path to the input FASTQ file containing sequencing reads.
    :param input_fasta: Path to the input FASTA file containing sequencing reads.
    :param output_tsv: Path to the output TSV file to store extracted UMIs. Defaults to the input filename with `.tsv` extension.
    :param umi_5_pattern: Regex pattern for the 5'-end UMI sequence. Defaults to a predefined pattern.
    :param umi_3_pattern: Regex pattern for the 3'-end UMI sequence. Defaults to a predefined pattern.
    :param n_errors: Number of errors allowed when searching for UMIs using fuzzy matching. Defaults to 2.
    :param not_found: Placeholder value to use when a UMI is not found in a read. Defaults to 'N/A'.
    :param unique_by_revcomp: If True, ensures unique UMI pairs by lexicographic ordering with respect to their reverse complement. Defaults to True.
    """
    
    # Initializing dataframe
    df = pd.DataFrame(columns=["5'-UMI", "3'-UMI"],
                      dtype=str)
    
    # Parsing input file
    if input_fastq is not None:
        reads = SeqIO.parse(input_fastq, 'fastq')
        if output_tsv is None:
            output_tsv = os.path.basename(input_fastq).replace('.fastq', '.tsv')
    elif input_fasta is not None:
        reads = SeqIO.parse(input_fasta, 'fasta')
        if output_tsv is None:
            output_tsv = os.path.basename(input_fasta).replace('.fasta', '.tsv')
    else:
        raise ValueError("Must pass either `input_fastq` or `input_fasta`!")
    
    # Parsing FASTQ file
    for read in reads:
        
        # Extracting UMIs
        umi_5, umi_3 = umis(str(read.seq), umi_5_pattern, umi_3_pattern, n_errors)
        df.loc[read.id] = [umi_5, umi_3]
    
    # Mapping to unique UMI pair by lexicographic priority with respect to reverse complement
    if unique_by_revcomp:
        df.fillna('', inplace=True)
        umis_for = df.values
        umis_rev = np.flip(df.applymap(lambda seq: str(Seq(seq).reverse_complement())).values, axis=1)
        for i in range(len(df)):
            df.iloc[i] = umis_for[i] if umis_for[i][0] < umis_rev[i][0] else umis_rev[i]
    
    # Setting placeholder value for not-found UMIs
    df.fillna(not_found, inplace=True)
    
    # Writing to output TSV
    if output_tsv is None:
        output_tsv = '.'.join(input_fa.split('.')[:-1]) + '.tsv'
    df.reset_index().to_csv(output_tsv, sep='\t', index=False)