#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import pandas as pd
from typing import Literal
from Bio import SeqIO


# ················································································· #

def fasta_to_df(input_fasta: str,
                field_delim: str='=',
                none_value: str='N/A') -> pd.DataFrame:
    """
    Parses a FASTA file into a ``pandas.DataFrame`` object. Makes header fields into columns.

    FASTA header must be formatted with space-delimited fields and a common name-value delimiter (default ``=``) for each field:

    ``>[ID] [field1]=[value1] [field2]=[value2] [field3]=[value3]``

    :param input_fasta: Path for the FASTA file to parse
    :param field_delim: The delimiter between a field name and value
    :param none_value: The empty value for fields (will be replaced with ``None``)
    :return: A DataFrame of parsed FASTA records.
    """

    # Looping over FASTA records
    tab_records = []
    fasta = SeqIO.parse(input_fasta, format='fasta')
    for record in fasta:

        # Parsing header
        header = record.description.split(' ')
        id = header.pop(0)
        fields = {'sequence': str(record.seq)}
        fields |= {name: value for name, value in [field.split(field_delim) for field in header]}

        # Formatting as tabular record
        record = pd.Series(name=id, data=fields)
        record.replace(none_value, None, inplace=True)
        tab_records.append(record)

    # Returning as DataFrame
    return pd.DataFrame(tab_records)


def fasta_to_tabular(output_file: str=None,
                     file_delim: Literal['tsv', 'csv']='csv',
                     **kwargs) -> None:
    """
    Takes a FASTA file, reformats it into a tabular delimited file.
    
    Makes header fields into columns.

    :param output_file: Path to output file. If `None`, reuses basename of ``input_fasta``.
    :param file_delim: The delimiter format to use for tabular format.
    :param kwargs: Passed to ``fasta_to_df()`` - see function for docs.
    """

    # Parse FASTA into a DataFrame
    df = fasta_to_df(**kwargs)

    # Reusing input file name, if no --output-file is passed
    if output_file is None:
        output_file = kwargs['input_fasta'].replace('fasta', file_delim)

    save_df(df, output_file, file_delim)


def save_df(df,
            output_file: str = None,
            file_delim: Literal['tsv', 'csv']='csv') -> None:
    """
    Takes a pandas DataFrame, saves it to either a 'csv' or 'tsv' file.

    :param df: DataFrame to save.
    :param output_file: Path to output file.
    :param file_delim: The delimiter format to use for tabular format.
    """

    # Saving to delimited file
    if file_delim == 'tsv':
        df.to_csv(output_file, sep='\t', index_label='id')
    if file_delim == 'csv':
        df.to_csv(output_file, sep=',', index_label='id')


# ················································································· #

if __name__ == '__main__':

    # Parsing input arguments
    parser = argparse.ArgumentParser(
        prog='translate-sam',
        usage='''
            python fasta-to-tabular.py
                [-h]
                -i/--input-fasta        INPUT.SAM [INPUT2.SAM ...]
                -o/--output-file        OUTPUT.[CSV/TSV/...]
                [-d/--file-delim        ['tsv', 'csv']           ]
                [-f/--field-delim       ['=']                    ]
                [-x/--none-value        ['N/A']                  ]  
                ''',
        description="Parse a FASTA into a delimited tabular format.",
        epilog='// Frederik Espersen Knudsen, 2024')

    parser.add_argument('-i', '--input-fasta',
                        type=str,
                        required=True,
                        help="FASTA file to parse into delimited format.")
    parser.add_argument('-o', '--output-file',
                        type=str,
                        required=False,
                        default=None,
                        help="Path to output the tabular file to.")
    parser.add_argument('-d', '--file-delim',
                        type=str,
                        required=False,
                        default='tsv',
                        help=r"The delimiter format of the tabular file format. OPTIONS: [tsv, csv].")
    parser.add_argument('-f', '--field-delim',
                        type=str,
                        required=False,
                        default='=',
                        help="The delimiter between name and value in FASTA header fields.")
    parser.add_argument('-x', '--none-value',
                        type=str,
                        required=False,
                        default='N/A',
                        help="How None-values are represented in FASTA header fields.")

    args = parser.parse_args()

    # ············································································· #

    # Parsing FASTA to tabular
    fasta_to_tabular(**vars(args))