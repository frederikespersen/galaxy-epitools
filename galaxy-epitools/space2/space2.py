#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import os
import glob
import SPACE2


# ················································································· #

def space2_binning(input_pdbs: list[str],
                   output_tsv: str):
    
    # Running SPACE2
    clustered_dataframe = SPACE2.agglomerative_clustering(input_pdbs, cutoff=1.25, n_jobs=-1)

    # Saving results
    clustered_dataframe.to_csv(output_tsv, sep='\t')
    
# ················································································· #

if __name__ == "__main__":
    
    # Setting up parser
    parser = argparse.ArgumentParser(description="Cluster antibody structures by epitope with SPACE2.")
    parser.add_argument("--input-pdbs", type=str, nargs='+', required=True, help="A directory containing PDB structures modelled with ABodyBuilder2" )
    parser.add_argument("--output-tsv", type=str, default='space2.tsv', help="A filename to write the output TSV to.")

    # ············································································· #

    # Calculating abundance
    args = parser.parse_args()
    space2_binning(**vars(args))
