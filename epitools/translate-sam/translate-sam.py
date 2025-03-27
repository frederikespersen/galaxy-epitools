#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import os
import regex as re
from typing import Literal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pysam import AlignmentFile
from pysam.libcalignedsegment import AlignedSegment


# ················································································· #

def translate_aligned_read(read: AlignedSegment,
                           within_mapping_only: bool = False,
                           truncate: bool = False) -> str:
    """
    Translates an aligned SAM read according to the reading frame of the alignment.
    Returns a trailing-end translation, from the alignment start to the end of the
    read (i.e. the translation does not end, when the alignment ends).
    
    An optional parameter ``within_mapping_only`` specifies whether to only return the AA sequence contained within the mapping (i.e. not to continue translation after the mapped part has ended. NOTE: Mapping of antibodies to a reference sequence may end prematurely at CDRH3).
    
    An optional parameter ``truncate`` specifies whether to truncate AA sequences at first stop codon '*'.

    :param read: An aligned SAM read parsed with ``pysam.AlignmentFile()``.
    :param within_mapping_only: Whether to limit translation to within the mapped stretch of the reads.
    :param truncate: Whether to limit translation to the first stop codon '*'.
    :return: A translation of the read in the alignment reading frame.
    """
    # Retrieving sequence
    read_sequence = read.query_sequence # If the read mapped to the reverse template, query_sequence is already reverse-complemented (Line 1949 in pysam/libcalignedsegment.pyx)
    
    # Determining open reading frame
    orf_start = read.query_alignment_start
    if within_mapping_only:
        orf_end = read.query_alignment_end
        codon_trim = len(read_sequence[orf_start:orf_end]) % 3
        orf_sequence = read_sequence[orf_start:orf_end] if codon_trim == 0 else read_sequence[orf_start:orf_end-codon_trim]
    else:
        codon_trim = len(read_sequence[orf_start:]) % 3
        orf_sequence = read_sequence[orf_start:] if codon_trim == 0 else read_sequence[orf_start:-codon_trim]

    # Translating
    read_translation = str(Seq(orf_sequence).translate())
    
    # Optional truncation
    if truncate:
        read_translation = read_translation.split('*')[0]
    
    return read_translation


def translate_sam(input_sam: str,
                  output_fasta: str,
                  within_mapping_only: bool,
                  truncate: bool,
                  val_min_mapq: int=60) -> None:
    """
    Takes a SAM of reads aligned to a template, writes a FASTA of reads translated
    within the open reading frame of their alignment.

    The translations are validated based on:
        1) Mapping quality.
    
    An optional parameter ``within_mapping_only`` specifies whether to only return the AA sequence contained within the mapping (i.e. not to continue translation after the mapped part has ended. NOTE: Mapping of antibodies to a reference sequence may end prematurely at CDRH3).
    
    An optional parameter ``truncate`` specifies whether to truncate AA sequences at first stop codon '*'.

    FASTA entries are returned with an ID corresponding to the original read ID, and a description denoting which reference the read was mapped to ``mapped_to``.

    :param input_sam: Path for the input SAM.
    :param output_fasta: Path for the output FASTA. If `None`, reuses basename of ``input_fasta``  (i.e. ``sample1.sam`` => ``sample1.fasta``).
    :param within_mapping_only: Whether to limit translation to within the mapped stretch of the reads.
    :param truncate: Whether to limit translation to the first stop codon '*'.
    :param val_min_mapq: The minimum ``minimap2`` mapping quality to accept (Maximum ``60``).
    """

    # Initializing a FASTA record container
    records = []

    # Looping over aligned reads
    with AlignmentFile(input_sam, 'r') as sam:
        for read in sam:

            # Skipping badly mapped reads
            if not read.is_mapped or read.mapping_quality < val_min_mapq:
                continue

            # Skipping reads that do not match template FRW4
            aa_seq = translate_aligned_read(read, within_mapping_only, truncate)

            # Formatting translation as FASTA record
            record = SeqRecord(Seq(aa_seq),
                               id=read.query_name,
                               description=f"mapped_to={read.reference_name}")
            records.append(record)
            
    # Reusing input file name, if no --output-fasta is passed
    if output_fasta is None:
        output_fasta = os.path.basename(input_sam).replace('.sam', '.fasta')
        
    # Saving records to FASTA
    if records:
        with open(output_fasta, 'w') as fasta:
            SeqIO.write(records, fasta, "fasta")


# ················································································· #

if __name__ == '__main__':

    # Getting commandline arguments
    parser = argparse.ArgumentParser(
        prog='translate-sam',
        usage='''
        python translate-sam.py
            [-h]
            -i/--input-sam          INPUT.SAM
            [-o/--output-fasta      OUTPUT.FASTA             ]
            [--truncate                                      ]   
            [--within-mapping-only                           ]   
            [-m/--val-min-mapq      [60]                     ]   
            ''',
        description="Translate single-chain Ab reads aligned to templates.",
        epilog='// Frederik Espersen Knudsen, 2024')

    parser.add_argument('-i', '--input-sam',
                        type=str,
                        required=True,
                        help="SAM file with reads aligned to templates.")
    parser.add_argument('-o', '--output-fasta',
                        type=str,
                        required=False,
                        default=None,
                        help="Path to output the translated FASTA to. Defaults to input path basename.")
    parser.add_argument('--truncate',
                        action='store_true',
                        required=False,
                        help="Truncate translation at the first stop codon '*'.")
    parser.add_argument('--within-mapping-only',
                        action='store_true',
                        required=False,
                        help="Limit translation to within the mapped stretch of the reads.")
    parser.add_argument('-m', '--val-min-mapq',
                        type=int,
                        required=False,
                        default=60,
                        help="The minimum minimap2 mapping quality to accept (Maximum 60).")

    args = parser.parse_args()

    # ············································································· #

    # Translating SAMs sequentially
    translate_sam(**vars(args))
