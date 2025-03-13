#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import os
import pandas as pd
from typing import Literal
from ImmuneBuilder import ABodyBuilder2, NanoBodyBuilder2


# ················································································· #

def predict(model: Literal["ABodyBuilder2", "NanoBodyBuilder2"],
            input_hc: str,
            input_lc: str = None,
            output_path: str = 'model.pdb'):
    
    # Specifying predictor model
    if model == "NanoBodyBuilder2":
        predictor = NanoBodyBuilder2()
        chains = {'H': input_hc}
        structure = predictor.predict(chains)
    elif model == "ABodyBuilder2":
        predictor = ABodyBuilder2()
        chains = {'H': input_hc, 'L': input_lc}
        structure = predictor.predict(chains)
    
    # Saving to PDB
    structure.save(output_path)

# ················································································· #

def predict_single(model: Literal["ABodyBuilder2", "NanoBodyBuilder2"],
                   input_hc: str,
                   input_lc: str = None,
                   output_dir: str = 'structures'):
    
    # Setting arbitrary PDB name
    output_path = os.path.join(output_dir, 'model.pdb')
    
    # Predicting and saving structure
    predict(model, input_hc, input_lc, output_path)

# ················································································· #

def predict_multiple(model: Literal["ABodyBuilder2", "NanoBodyBuilder2"],
                     input_tsv: str,
                     id_col: str,
                     hc_col: str,
                     lc_col: str,
                     output_dir: str):
    
    # Loading data
    df = pd.read_csv(input_tsv, sep='\t')
    
    # Extracting columns
    ids = df[id_col]
    input_hcs = df[hc_col]
    if model == "ABodyBuilder2":
        input_lcs = df[lc_col]
    elif model == "NanoBodyBuilder2":
        input_lcs = [None] * len(input_hcs)
    
    # Looping over entries
    for id, input_hc, input_lc in zip(ids, input_hcs, input_lcs):
        
        # Setting PDB name as ID
        output_path = os.path.join(output_dir, f'{id}.pdb')
        
        # Predicting and saving structure
        predict(model, input_hc, input_lc, output_path)

# ················································································· #

if __name__ == "__main__":

    # Setting up parsers
    parser = argparse.ArgumentParser(description="Predict antibody structures using ImmuneBuilder.")
    subparsers = parser.add_subparsers(dest="mode", required=True)
    
    single_parser = subparsers.add_parser("single", help="Predict a single antibody structure")
    single_parser.add_argument("--model", choices=["ABodyBuilder2", "NanoBodyBuilder2"], required=True)
    single_parser.add_argument("--input-hc", required=True, help="Heavy chain sequence")
    single_parser.add_argument("--input-lc", help="Light chain sequence (only for ABodyBuilder2)")
    single_parser.add_argument("--output-dir", required=True, help="Output directory")
    
    multiple_parser = subparsers.add_parser("multiple", help="Predict multiple antibody structures")
    multiple_parser.add_argument("--model", choices=["ABodyBuilder2", "NanoBodyBuilder2"], required=True)
    multiple_parser.add_argument("--input-tsv", required=True, help="TSV file with antibody sequences")
    multiple_parser.add_argument("--id-col", type=str, required=False, default='id', help="Name of the ID column [Default 'ID']")
    multiple_parser.add_argument("--hc-col", type=str, required=False, default='hc', help="Name of the Heavy Chain column [Default 'HC']")
    multiple_parser.add_argument("--lc-col", type=str, required=False, default='lc', help="Name of the Light Chain column [Default 'LC']")
    multiple_parser.add_argument("--output-dir", required=True, help="Output directory")
    
    # ················································································· #
    
    # Parsing arguments and creating output directory
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Predicting either a single or multiple structures
    if args.mode == "single":
        del args.mode
        predict_single(**vars(args))
    elif args.mode == "multiple":
        del args.mode
        predict_multiple(**vars(args))

