#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import pandas as pd
from typing import Literal


# ················································································· #

def extract_ab_region(numbering_file: str,
                      name: str,
                      start: str,
                      end: str,
                      output_file: str=None,
                      gap_symbol: str='-') -> None:

    # Loading numbering file
    numbering = pd.read_table(numbering_file, index_col=0, sep=None, engine='python')

    # Setting index as FASTA IDs (instead of entire header)
    fasta_id = lambda header: header.split(' ')[0]
    numbering.index = numbering.index.map(fasta_id)
    numbering.index.name = 'id'

    # Subsetting region
    start_idx = numbering.columns.get_loc(start)
    end_idx = numbering.columns.get_loc(end)
    region = numbering.iloc[:, start_idx:end_idx+1].sum(axis=1).str.replace(gap_symbol, '')
    region.name = name

    # Saving
    if output_file is None:
        output_file = f'{name}.tsv'
    file_delim = output_file.split('.')[-1]
    if file_delim not in ['tsv', 'csv']:
        raise ValueError(f"File format must be either '.tsv' or '.csv'! (.{file_delim} was passed for file {output_file})")
    save_df(region.to_frame(),
            output_file,
            file_delim)


def save_df(df,
            output_file: str,
            file_delim: Literal['tsv', 'csv'] = 'csv') -> None:
    """
    Takes a pandas DataFrame, saves it to either a 'csv' or 'tsv' file.

    :param df: DataFrame to save.
    :param output_file: Path to output file.
    :param file_delim: The delimiter format to use for tabular format.
    """

    # Saving to delimited file
    if file_delim == 'tsv':
        df.to_csv(output_file, sep='\t', index_label='id')
    elif file_delim == 'csv':
        df.to_csv(output_file, sep=',', index_label='id')


# ················································································· #

if __name__ == '__main__':

    # Parsing input arguments
    parser = argparse.ArgumentParser(
        prog='extract-ab-region',
        usage='''
            python extract-ab-region.py
                [-h]
                -i/--numbering-file     NUMBERING.[CSV]
                -n/--name               NAME
                -s/--start              START
                -e/--end                END
            [   -o/--output-file        OUTPUT.[TSV]    ]
            [   -g/--gap-symbol         [-]             ]
                ''',
        description="Extract Ab region sequences from a CSV of numbered positions (the output of ``parse-anarci.py``).",
        epilog='// Frederik Espersen Knudsen, 2024')

    parser.add_argument('-i', '--numbering-file',
                        type=str,
                        required=True,
                        help="Path to CSV-file containing columns of numbered positions. First column will be used as index.")
    parser.add_argument('-n', '--name',
                        type=str,
                        required=False,
                        default='Region',
                        help="The name to use for the region in the output file.")
    parser.add_argument('-s', '--start',
                        type=str,
                        required=True,
                        help="The numbered position where the region starts (inclusive).")
    parser.add_argument('-e', '--end',
                        type=str,
                        required=True,
                        help="The numbered position where the region ends (inclusive).")
    parser.add_argument('-o', '--output-file',
                        type=str,
                        required=False,
                        default=None,
                        help="The path to write the output file to. (Defaults to [region].tsv)")
    parser.add_argument('-g', '--gap-symbol',
                        type=str,
                        required=False,
                        default='-',
                        help="The gap symbol used in the input numbering file (ANARCI default is '-').")

    args = parser.parse_args()

    # ············································································· #

    # Parsing ANARCI output files
    extract_ab_region(**vars(args))