#!/usr/bin/env python3

# ················································································· #

# Imports
import argparse
from collections.abc import Callable
from typing import Literal
import pandas as pd
import numpy as np
import Levenshtein


# ················································································· #

MEASURES = {
    "levenshtein": Levenshtein.distance,
}

REPRESENTATIONS = {
    'standard': {},
    'specifica': {
        'A': 'C', # Hydrophobic, Small
        'C': 'H', # Cysteine
        'D': 'A', # Acidic
        'E': 'A', # Acidic
        'F': 'D', # Aromatic
        'G': 'J', # Glycine
        'H': 'F', # Basic
        'I': 'B', # Hydrophobic, Aliphatic
        'K': 'F', # Basic
        'L': 'B', # Hydrophobic, Aliphatic
        'M': 'K', # Methionine
        'N': 'E', # Amide
        'P': 'I', # Proline
        'Q': 'E', # Amide
        'R': 'F', # Basic
        'S': 'G', # Hydroxyl
        'T': 'G', # Hydroxyl
        'V': 'B', # Hydrophobic, Aliphatic
        'W': 'D', # Aromatic
        'Y': 'D', # Aromatic
    }
}


# ················································································· #

def calculate_sequence_distance_matrix(query_seqs: np.ndarray[str],
                                       reference_seqs: np.ndarray[str],
                                       measure: Callable):
    
    # Check equal amounts of sequence subsets
    if len(query_seqs.shape) == 2:
        assert query_seqs.shape[1] == reference_seqs.shape[1], "Input arrays must have same second dimension, i.e. number of subsequences."
    assert len(query_seqs.shape) <= 2, "Input sequence arrays must be of dimension [n,m], where n is number of sequence sets and m number of subsequences in a set."
    
    # Find unique sequences and their indices
    unique_query_seqs, query_idx = np.unique(query_seqs, return_inverse=True)
    unique_reference_seqs, reference_idx = np.unique(reference_seqs, return_inverse=True)
    
    # Initialize unique distance matrix
    q = len(unique_query_seqs)
    r = len(unique_reference_seqs)
    unique_distance_matrix = np.zeros((q, r), dtype=int)
    
    # Looping over unique pairs of sequence comparisons
    for i, q_seq in enumerate(unique_query_seqs):
        for j, r_seq in enumerate(unique_reference_seqs):
            unique_distance_matrix[i, j] = measure(q_seq, r_seq)
    
    # Filling distance matrix by mapping back pair distances by indices
    return unique_distance_matrix[np.reshape(query_idx, query_seqs.shape)[:,None],
                                  np.reshape(reference_idx[None,:], reference_seqs.shape)]

    
def sequence_distance_matrix(input_set1_tsv: str,
                             input_set2_tsv: str,
                             measure: Literal['levenshtein'],
                             representation: Literal['standard', 'specifica'],
                             output_tsv: str,
                             sum_subsets: bool=True,
                             input_set1_tsv_idcol: bool=False,
                             input_set2_tsv_idcol: bool=False):
    
    # Loading data
    input_set1 = pd.read_csv(input_set1_tsv, sep='\t', index_col=0 if input_set1_tsv_idcol else None).astype(str)
    input_set2 = pd.read_csv(input_set2_tsv, sep='\t', index_col=0 if input_set2_tsv_idcol else None).astype(str)
    
    # QC of data input
    assert (input_set1.columns == input_set2.columns).all(), "The input TSVs must have the same columns! (Case-sensitive headers)"
    
    # Transforming sequence tokenization
    representation = str.maketrans(REPRESENTATIONS[representation])
    input_set1 = input_set1.map(lambda seq: seq.translate(representation))
    input_set2 = input_set2.map(lambda seq: seq.translate(representation))
    set1, set2 = input_set1.values, input_set2.values
    
    # Calculating sequence distance matrix
    measure = MEASURES[measure]
    dists = calculate_sequence_distance_matrix(set1, set2, measure)
    
    # Formatting output
    if sum_subsets:
        dists = dists.sum(axis=2)
    else:
        dists = dists.tolist()
    output_dists = pd.DataFrame(dists,
                                columns=input_set2.index,
                                index=input_set1.index).T
    output_dists.index.rename(None, inplace=True)
    
    # Writing to TSV
    output_dists.to_csv(output_tsv, sep='\t')
    

# ················································································· #

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        prog='sequence-distance',
        description="""
        Compute the sequence distance matrix between two sets of sequences using a specified
        distance measure and representation mapping.
        """,
        epilog='// Generated Script'
    )
    
    parser.add_argument('--input-set1-tsv',
                        type=str,
                        required=True,
                        help="Path to the first input TSV file containing sequences.")
    parser.add_argument('--input-set2-tsv',
                        type=str,
                        required=True,
                        help="Path to the second input TSV file containing sequences.")
    parser.add_argument('--measure',
                        type=str,
                        choices=MEASURES.keys(),
                        required=True,
                        help="Distance measure to use (e.g., 'levenshtein').")
    parser.add_argument('--representation',
                        type=str,
                        choices=REPRESENTATIONS.keys(),
                        default='standard',
                        help="Sequence representation type (default: 'standard').")
    parser.add_argument('--output-tsv',
                        type=str,
                        required=True,
                        help="Path to save the computed sequence distance matrix as a TSV file.")
    parser.add_argument('--sum-subsets',
                        action='store_true',
                        help="Sum the distances over sequence subsets (default: False).")
    parser.add_argument('--input-set1-tsv-idcol',
                        action='store_true',
                        help="Treat the first column of input-set1 as an index column.")
    parser.add_argument('--input-set2-tsv-idcol',
                        action='store_true',
                        help="Treat the first column of input-set2 as an index column.")
    
    # ················································································· #

    args = parser.parse_args()
    
    sequence_distance_matrix(**vars(args))