#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import re
import pandas as pd
from typing import Literal
from Bio import SeqIO


# ················································································· #

def parse_liabilities(input_strs: list[str]) -> dict[list]:
    
    # Initializing dictionary for liability patterns
    liability_patterns = {}
    
    # Looping over input strings
    for input_str in input_strs:
        
        # Separating liability name and regex patterns, and separating different patterns
        liability, patterns = input_str.split('=')
        patterns = patterns.split(',')
        liability = liability.replace("'","")
        
        liability_patterns[liability] = patterns
        
    return liability_patterns


def find_liabilities(sequence: str,
                     liability_patterns: dict[list]):
    
    # Looping over liabilities
    liabilities = []
    for liability, patterns in liability_patterns.items():
        
        # Looping over patterns
        for pattern in patterns:
            
            # Finding matches to patterns
            for hit in re.finditer(pattern, sequence):
                                
                # Parsing match
                liabilities.append({
                    'liability': liability,
                    'match': hit.group(),
                    'start': hit.span()[0] + 1,    # Non-zero indexing for human interpretability
                    'end': hit.span()[1]
                })
                
    if len(liabilities) == 0:
        return pd.DataFrame(columns=['liability', 'match', 'start', 'end'])
    else:
        return pd.DataFrame(liabilities).sort_values(['start', 'end'], ignore_index=True)
    

def sequence_liabilities(id: str,
                         sequence: str,
                         liability_patterns: list[str]=[],
                         format: Literal['condensed', 'records', 'pretty'] = 'condensed',
                         padding: str='.',):
    
    # Removing newlines from sequence
    sequence = sequence.replace('\n','')
    
    # Parsing input strings for liabilities
    liability_patterns = parse_liabilities(liability_patterns)
    
    # Finding liabilities
    liabilities = find_liabilities(sequence, liability_patterns)
    
    # The standard format is single liabilities as records in a table
    if format == 'records':
        
        # Adding id column
        liabilities = pd.concat([
            pd.DataFrame({'id': [id]*len(liabilities)}),
            liabilities
        ], axis=1)
            
    # Condensed is a single-line in a table
    elif format == 'condensed':
        
        # Formatting liabilities as a list of tuples
        liabilities = [tuple(liabilities.loc[idx].values) for idx in liabilities.index]
        
        # Adding id column
        liabilities = pd.DataFrame({'id': [id], 'liabilities': [liabilities]})
            
    # Pretty displays the liabilities in an alignment
    elif format == 'pretty':
        
        # Making liabilities into position-wise alignments by padding
        align_liabilities = []
        for idx in liabilities.index:
            match = liabilities.loc[idx, 'match']
            pad_right = liabilities.loc[idx, 'start'] - 1
            pad_left = len(sequence) - liabilities.loc[idx, 'end']
            align_match = pad_right * padding + match + pad_left * padding
            align_liabilities.append(align_match)
        
        alignment = pd.concat([
            pd.DataFrame({'rowheader': [id], 'sequence': [sequence]}),
            pd.DataFrame({'rowheader': '  * ' + liabilities['liability'], 'sequence': align_liabilities})
        ])
        rowheader_len = alignment['rowheader'].map(len).max() + 3
        alignment['rowheader'] = alignment['rowheader'].map(lambda s: f"{s: <{rowheader_len}}")
        liabilities = alignment.sum(axis=1)
    
    return liabilities


# ················································································· #

if __name__ == '__main__':

    # Setting up parsers
    parser = argparse.ArgumentParser(description="Find sequence liabilities.")
    subparsers = parser.add_subparsers(dest="mode", required=True)
    
    single_parser = subparsers.add_parser("single", help="Find liabilities in a single sequence.")
    single_parser.add_argument("--input-sequence", type=str, required=True, help="An amino acid sequence to look for liabilities in.")
    single_parser.add_argument("--liability-patterns", nargs='+', type=str, help='Liabilities and a comma-separated list of regex patterns to discover them, provided in the format "'+'"'+"'[LiabilityName]':[RegexPattern1],[RegexPattern2]" + '"' + ". For instance "  + '"' + "'Asn Deamidation'=N[GNST],GN[FGTY]" + '"' + ".")
    single_parser.add_argument("--format", type=str, default='records', required=False, help="What format to return results in. Options are ['condensed', 'records', 'pretty']")
    single_parser.add_argument("--output-tsv", type=str, required=True, help="Path to TSV to write output to.")

    multiple_parser = subparsers.add_parser("multiple", help="Find liabilities in a set of sequences.")
    multiple_parser.add_argument("--input-fasta", type=str, required=True, help="A FASTA file with amino acid sequences to look for liabilities in.")
    multiple_parser.add_argument("--liability-patterns", nargs='+', type=str, help='Liabilities and a comma-separated list of regex patterns to discover them, provided in the format "'+'"'+"'[LiabilityName]':[RegexPattern1],[RegexPattern2]" + '"' + ". For instance "  + '"' + "'Asn Deamidation'=N[GNST],GN[FGTY]" + '"' + ".")
    multiple_parser.add_argument("--format", type=str, default='records', required=False, help="What format to return results in. Options are ['condensed', 'records', 'pretty']")
    multiple_parser.add_argument("--output-tsv", type=str, required=True, help="Path to TSV to write output to.")

    # ················································································· #
    
    # Parsing arguments and creating output directory
    args = parser.parse_args()
    
    # Finding liabilities in a single sequence and writing output to TSV
    if args.mode == "single":
        liabilities = sequence_liabilities(id='Input Sequence',
                                           sequence=args.input_sequence,
                                           liability_patterns=args.liability_patterns,
                                           format=args.format)
        liabilities.to_csv(args.output_tsv, sep='\t', header=False, index=False)
        
    # Finding liabilities in multiple sequences and writing output to TSV
    elif args.mode == "multiple":
        
        # Parsing sequences
        sequences = {}
        for record in SeqIO.parse(args.input_fasta, 'fasta'):
            sequences[record.id] = str(record.seq)
        
        # Determining liabilities for all sequences
        liability_dfs = []
        for id, sequence in sequences.items():
            liability_dfs.append(sequence_liabilities(id=id,
                                                      sequence=sequence,
                                                      liability_patterns=args.liability_patterns,
                                                      format=args.format))
        
        # Combining results and writing to output
        liabilities = pd.concat(liability_dfs)
        liabilities.to_csv(args.output_tsv, sep='\t', header=False, index=False)
