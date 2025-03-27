#!usr/bin/env python3

# ················································································· #

# Imports
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import Levenshtein
import argparse

# ················································································· #

def trim_frame_to_codons(sequence, start_idx):
    
    trimmed_length = len(sequence[start_idx:]) - (len(sequence[start_idx:]) % 3)
    return sequence[start_idx:start_idx + trimmed_length]


def translate_orfs(sequence: Seq) -> pd.Series:
    
    # Computing reverse complement
    sequence_for = sequence
    sequence_rev = sequence.reverse_complement()
    
    # Translating all ORFs end-to-end
    orfs = pd.Series({
        "orf_forward_1": trim_frame_to_codons(sequence_for, 0).translate(),
        "orf_forward_2": trim_frame_to_codons(sequence_for, 1).translate(),
        "orf_forward_3": trim_frame_to_codons(sequence_for, 2).translate(),
        "orf_reverse_1": trim_frame_to_codons(sequence_rev, 0).translate(),
        "orf_reverse_2": trim_frame_to_codons(sequence_rev, 1).translate(),
        "orf_reverse_3": trim_frame_to_codons(sequence_rev, 2).translate()
    }).astype(str)
    
    return orfs


def load_fasta(fasta_path: str) -> pd.Series:
    
    entries = {}
    for entry in SeqIO.parse(fasta_path, 'fasta'):
        entries[entry.id] = str(entry.seq)
    return pd.Series(entries).rename("id")


def best_template_match(orfs: pd.DataFrame,
                        templates: pd.Series) -> pd.DataFrame:
    
    # Determining distance matrix between templates and orfs
    distances = []
    for template_sequence in templates:
        distances.append(orfs.map(lambda orf: Levenshtein.distance(orf, template_sequence)).values)
    distances = np.stack(distances)

    # Determining best template for each read
    template_match_idx = np.argmin(distances.min(axis=2), axis=0)

    # Determining best ORF for each read
    orf_match_idx = np.argmin(distances[template_match_idx, np.arange(len(orfs)), :], axis=1)

    # Returning best ORFs and best templates
    orf_match = pd.Series(orfs.values[np.arange(len(orfs)), orf_match_idx], dtype=str, index=orfs.index)
    template_match = pd.Series(templates.index.values[template_match_idx], dtype=str, index=orfs.index)
    
    # Assembling dataframe
    matches = pd.DataFrame({'best_match_translation': orf_match,
                            'best_match_template_name': template_match})
    return matches
    

def translate_orfs_fa(fa_path: str) -> pd.DataFrame:
    
     # Initiliazing list
    orfs = []
    
    # Initializing file parser
    try:
        entries = SeqIO.parse(fa_path, 'fastq')
    except: 
        try:
            entries = SeqIO.parse(fa_path, 'fasta')
        except:
            raise ValueError("Must provide a file in either FASTA or FASTQ format!")
    
    # Translating ORFs for each read
    for seqrecord in entries:
        orfs.append(translate_orfs(seqrecord.seq).rename(seqrecord.id))
        
    # Assembling dataframe
    orfs = pd.DataFrame(orfs)
    orfs.index.name = "id"
    return orfs

# ················································································· #

if __name__ == "__main__":
    
    # Setting up parser
    parser = argparse.ArgumentParser(description="Translate all open reading frames.")
    parser.add_argument("--input-fa", type=str, required=True, help="A FASTA or FASTQ file of DNA sequences to translate" )
    parser.add_argument("--templates-fasta", type=str, required=False, help="A FASTA of amino acid sequence templates. Used to select the best ORF for each read, by minimum edit distance. If multiple templates are provided, the template-ORF pair with the lowest edit distance is returned.")
    parser.add_argument("--output-tsv", type=str, required=False, default='translations.csv', help="A filename to write the output TSV to.")

    # ············································································· #
    
    args = parser.parse_args()

    # Translating ORFs
    orfs = translate_orfs_fa(fa_path=args.input_fa)
    
    # Finding best match if templates are available
    if args.templates_fasta is not None:    
        templates = load_fasta(args.templates_fasta)
        matches = best_template_match(orfs, templates)
        output = pd.concat([orfs, matches], axis=1)
    else:
        output = orfs
    
    # Writing to TSV
    output.to_csv(args.output_tsv, sep='\t')