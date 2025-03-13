#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import pandas as pd
from collections.abc import Iterable
from typing import Literal
import itertools


# ················································································· #

def load_ngs_dataset(input_tsv: str,
                     framework_col: str,
                     cdr_cols: list[str],
                     umi_cols: list[str] = []) -> pd.DataFrame:
    
    # Loading data
    data = pd.read_csv(input_tsv, delimiter='\t')
    
    # Assembling data
    df = pd.concat([
        pd.DataFrame({'framework': data[framework_col]}),
        data[cdr_cols],
        data[umi_cols]
        ], axis=1)
    
    return df, (cdr_cols, umi_cols)

# ················································································· #

def abundance_ngs_dataset(output_tsv: str = None,
                          normalize: int = False,
                          **load_ngs_dataset_kwargs) -> pd.DataFrame:
    
    # Loading data
    df, (cdr_cols, umi_cols) = load_ngs_dataset(**load_ngs_dataset_kwargs)
    
    # Removing redundant UMI pairs (where avilable)
    if umi_cols != []:
        has_umis = df[umi_cols].notna().all(axis=1)
        df = pd.concat([
            df[~has_umis],
            df[has_umis].drop_duplicates(subset=umi_cols)
        ]).drop(umi_cols, axis=1)
    
    # Remove hits with empty CDRs
    df = df.dropna(subset=cdr_cols)
    
    # Counting reads per hits
    df = df.value_counts().reset_index().rename({'count': 'Count'}, axis=1)
    
    # Determining primary framework (assumed to be library) and removing contaminants
    framework_counts = df.value_counts('framework')
    primary_framework = framework_counts.index[framework_counts.argmax()]
    df = df[df['framework'] == primary_framework]
    
    # Calculating normalized counts
    if normalize:
        df['Count normalized'] = df['Count'] / normalize
    
    # Either returning as DataFrame object or saving to file
    if output_tsv is None:
        return df
    else:
        df.to_csv(output_tsv, sep='\t', index=False)

# ················································································· #

if __name__ == "__main__":
    
    # Setting up parser
    parser = argparse.ArgumentParser(description="Calculate antibody abundance in a NGS dataset.")
    parser.add_argument("--input-tsv", type=str, required=True, help="Path to the input TSV dataset")
    parser.add_argument("--framework-col", type=str, required=False, default='framework', help="Name of the framework column [Default 'framework']")
    parser.add_argument("--cdr-cols", type=str, nargs='+', required=False, default=['cdrl1','cdrl2','cdrl3','cdrh1','cdrh2','cdrh3'], help="Names of CDR columns [Default ['cdrl1','cdrl2','cdrl3','cdrh1','cdrh2','cdrh3']]")
    parser.add_argument("--umi-cols", type=str, nargs='+', required=False, default=[], help="Names of UMI columns [If not provided, no UMI-bias correction is performed]")
    parser.add_argument("--output-tsv", default='ab-abundance.tsv', help="Path to save the processed dataset to as TSV")
    
    # ············································································· #

    # Calculating abundance
    args = parser.parse_args()
    abundance_ngs_dataset(**vars(args))