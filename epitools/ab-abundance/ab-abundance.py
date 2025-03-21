#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import pandas as pd


# ················································································· #

def load_ngs_dataset(input_tsv: str,
                     framework_col: str,
                     cdr_cols: list[str],
                     umi_cols: list[str] = []) -> pd.DataFrame:
    
    # Loading data
    data = pd.read_csv(input_tsv, delimiter='\t', index_col=0)
    
    # Assembling data
    df = pd.concat([
        data[framework_col].rename('framework'),
        data[cdr_cols],
        data[umi_cols]
        ], axis=1)
    
    return df, (cdr_cols, umi_cols)

# ················································································· #

def abundance_ngs_dataset(output_tsv: str = None,
                          normalize: int = False,
                          decontaminate: bool = False,
                          return_ids: bool=False,
                          **load_ngs_dataset_kwargs) -> pd.DataFrame:
    
    # Loading data
    df, (cdr_cols, umi_cols) = load_ngs_dataset(**load_ngs_dataset_kwargs)
    
    # Removing redundant UMI pairs (where available)
    if umi_cols != []:
        has_umis = df[umi_cols].notna().all(axis=1)
        df = pd.concat([
            df[~has_umis],
            df[has_umis].drop_duplicates(subset=umi_cols)
        ]).drop(umi_cols, axis=1)
    
    # Remove hits with empty CDRs
    df = df.dropna(subset=cdr_cols)
    
    # Counting reads per hits
    df_counts = df.value_counts()

    # Optionally appending IDs
    if return_ids:
        df.index.rename('id', inplace=True)
        df_ids = df.reset_index().groupby([*df.columns])['id'].apply(list).rename('ids')
        df = pd.merge(df_counts, df_ids, left_index=True, right_index=True).reset_index()
    else: 
        df = df_counts.reset_index()
    
    # Determining primary framework (assumed to be library) and removing contaminants
    if decontaminate:
        framework_counts = df.value_counts('framework')
        primary_framework = framework_counts.index[framework_counts.argmax()]
        df = df[df['framework'] == primary_framework]
    
    # Calculating normalized counts
    if normalize:
        df['count_normalized'] = df['count'] / normalize
    
    # Either returning as DataFrame object or saving to file
    if output_tsv is None:
        return df
    else:
        df.to_csv(output_tsv, sep='\t', index=False)

# ················································································· #

if __name__ == "__main__":
    
    # Setting up parser
    parser = argparse.ArgumentParser(description="Calculate antibody abundance")
    parser.add_argument("--input-tsv", type=str, required=True, help="Path to the input TSV dataset")
    parser.add_argument("--framework-col", type=str, required=False, default='framework', help="Name of the framework column [Default 'framework']")
    parser.add_argument("--cdr-cols", type=str, nargs='+', required=False, default=['cdrl1','cdrl2','cdrl3','cdrh1','cdrh2','cdrh3'], help="Names of CDR columns [Default ['cdrl1','cdrl2','cdrl3','cdrh1','cdrh2','cdrh3']]")
    parser.add_argument("--umi-cols", type=str, nargs='+', required=False, default=[], help="Names of UMI columns [If not provided, no UMI-bias correction is performed]")
    parser.add_argument("--decontaminate", action='store_true', required=False, help="Whether to remove entires of all but the highest abundance framework, assuming they are contaminants.")
    parser.add_argument("--return-ids", action='store_true', required=False, help="Whether to return a column of lists of IDs of the entries that contributed to the count (For backtracking to wells).")
    parser.add_argument("--output-tsv", default='ab-abundance.tsv', help="Path to save the processed dataset to as TSV")
    
    # ············································································· #

    # Calculating abundance
    args = parser.parse_args()
    abundance_ngs_dataset(**vars(args))