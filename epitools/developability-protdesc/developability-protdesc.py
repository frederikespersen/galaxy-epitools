#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import pandas as pd


# ················································································· #

def calc_tab(schro_prot_desc: pd.DataFrame) -> pd.Series:
    """
    Takes a Protein Descriptors dataframe,
    returns the  Thera-SAbDaB (TAB) developability measure from the Charlotte Deane lab.

    :param: schro_prot_desc: A Protein Descriptors dataframe, generated by Schrodinger through the Calculate Protein Descriptors Task.
    :returns: A series of TAB values for each entry in the dataframe.
    """
    
    #TODO: Find TAB Model
    return 


def calc_tada(schro_prot_desc: pd.DataFrame) -> pd.Series:
    """
    Takes a Protein Descriptors dataframe,
    returns the Therapeutic Antibody Developability (TA-DA developability) measure from the XXX lab.

    :param: schro_prot_desc: A Protein Descriptors dataframe, generated by Schrodinger through the Calculate Protein Descriptors Task.
    :returns: A series of TAB values for each entry in the dataframe.
    """
    
    #TODO: Find TA-DA Model
    return 


def developability_measures(schroprotdesc_csv: str,
                            output_tsv: str=None) -> pd.DataFrame:
    """
    Takes the path to a Protein Descriptors CSV file,
    returns or writes the dataframe with selected descriptors and calculated developability measures as a TSV.

    :param: schroprotdesc_csv: Path to a Protein Descriptors CSV, generated by Schrodinger through the Calculate Protein Descriptors Task.
    :returns: If no ``output_tsv`` is specified, a curated dataframe of developability measures. If ``output_tsv`` is specified, this dataframe is written to file.
    """
    # Loading data
    schro_prot_desc = pd.read_csv(schroprotdesc_csv, index_col=0)
    
    # Curating select descriptors
    developability_prot_descs = ['All_AggScore']
    developability = pd.concat([
        schro_prot_desc[developability_prot_descs],
        calc_tab(schro_prot_desc).to_frame(),
        calc_tada(schro_prot_desc).to_frame()
    ])
    
    if output_tsv is None:
        return developability
    else:
        developability.to_csv(output_tsv, sep='\t')


# ················································································· #

if __name__ == '__main__':
    
    # Parsing input arguments
    parser = argparse.ArgumentParser(
        prog='extract-ab-region',
        usage='''
            python extract-ab-region.py
                [-h]
                -i/--schroprotdesc-csv  PROTEIN_DESCRIPTORS.CSV
            [   -o/--output-tsv        OUTPUT.TSV              ]
                ''',
        description="Calculate and curate developability measures from a Schrodinger Protein Descriptors file (See of ``schrodinger-prot-desc``).",
        epilog='// Frederik Espersen Knudsen, 2025')

    parser.add_argument('-i', '--schroprotdesc-csv',
                        type=str,
                        required=True,
                        help="Path to a Protein Descriptors dataframe, generated by Schrodinger through the Calculate Protein Descriptors Task.")
    parser.add_argument('-o', '--output-tsv',
                        type=str,
                        required=False,
                        default='developability.tsv',
                        help="The path to write the output TSV file to. (Defaults to developability.tsv)")

    args = parser.parse_args()

    # ············································································· #

    # Parsing ANARCI output files
    developability_measures(**vars(args))