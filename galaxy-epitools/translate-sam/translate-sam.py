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

def translate_aligned_read(read: AlignedSegment) -> str:
    """
    Translates an aligned SAM read according to the reading frame of the alignment.
    Returns a trailing-end translation, from the alignment start to the end of the
    read (i.e. the translation does not end, when the alignment ends).

    :param read: An aligned SAM read parsed with ``pysam.AlignmentFile()``.
    :return: A translation of the read in the alignment reading frame.
    """
    read_sequence = read.query_sequence
    orf_start = read.query_alignment_start
    codon_trim = len(read_sequence[orf_start:]) % 3
    orf_sequence = read_sequence[orf_start:] if codon_trim == 0 else read_sequence[orf_start:-codon_trim]
    return str(Seq(orf_sequence).translate())


def load_validation_end_patterns(fasta: str,
                                 n_last: int,
                                 n_sub: int,
                                 fasta_type: Literal['dna', 'aa']= 'dna') -> dict[str, str]:
    """
    Assembles substitution-sensitive regex patterns from the C-terminal-ends of
    loaded template sequences.

    :param fasta: Template sequences, from which to generate C-ter regex patterns.
    :param n_last: How many of the C-ter template residues to use for patterns.
    :param n_sub: How many substitutions to tolerate in the regex pattern.
    :param fasta_type: Whether the templates are written as nucleotides or amino acids.
    :return: A dictionary of a C-ter regex pattern per loaded template.
    """
    # Checking that any end residues are specified for validation
    assert n_last > 0, "Argument `n_last` must be ≥ 0!"

    # Loading templates, and translating if DNA
    fasta = SeqIO.parse(fasta, format='fasta')
    if fasta_type == 'dna':
        templates = {record.id: str(Seq(record.seq).translate()) for record in fasta}
    elif fasta_type == 'aa':
        templates = {record.id: str(record.seq) for record in fasta}
    else:
        raise ValueError("Specify whether FASTA templates are in nucleotides `dna` or amino acids `aa`!")

    # Extracting end residues
    validation_end_seqs = {name: seq[-n_last:] for name, seq in templates.items()}

    # Generating patterns
    validation_patterns = {name: end_pattern(seq, n_sub) for name, seq in validation_end_seqs.items()}
    return validation_patterns


def end_pattern(pattern: str,
                n_sub: int) -> str:
    """
    Takes a regex pattern, assembles a derived regex pattern that tolerates
    a set amount of substitutions, and which returns a group including the
    original pattern and any leading text in a matching string.

    If there is no match within the substitution threshold, this regex
    pattern return no groups (i.e. ``re.match(str,pattern)`` returns ``None``).

    :param pattern: A regex pattern.
    :param n_sub: How many substitutions to tolerate for match to ``pattern``.
    :return: A substitution-sensitive derived regex pattern that returns
    matches including leading strings.
    """
    return "(^.*" + pattern + ")"+"{s<=" + str(n_sub) + "}"


def umis(seq: str,
         umi_5_pattern: str=r'T{3}[ACG]{4}T{2}[ACG]{4}T{2}[ACG]{4}T{2}[ACG]{4}T{3}',
         umi_3_pattern: str=r'A{3}[TCG]{4}A{2}[TCG]{4}A{2}[TCG]{4}A{2}[TCG]{4}A{3}',
         n_errors: int=2,
         not_found=None) -> tuple:
    """
    Takes a sequence, tries to find UMIs within a certain margin of error,
    and returns those UMIs.

    :param seq: A sequence in which to search for the UMI in.
    :param umi_5_pattern: The pattern for a 5'-end UMI.
    :param umi_3_pattern: The pattern for a 3' UMI
    :param n_errors: How many errors to accept when searching for either UMI.
    :param not_found: What to return, when a UMI is not found.
    :return: A tuple of (5' UMI, 3' UMI). If UMIs are not found, value will be ``None``.
    """
    umi_5_pattern = '(' + umi_5_pattern + '){s<=' + str(n_errors) + '}'
    umi_3_pattern = '(' + umi_3_pattern + '){s<=' + str(n_errors) + '}'
    umi_matches = (re.search(pat, seq) for pat in [umi_5_pattern, umi_3_pattern])
    return tuple(match.group() if match is not None else not_found for match in umi_matches)


def translate_sam(input_sam: str,
                  output_fasta: str,
                  templates_fasta: str,
                  fasta_type: Literal['dna', 'aa']= 'dna',
                  val_min_mapq: int=60,
                  val_n_last: int=1,
                  val_n_sub: int=1) -> None:
    """
    Takes a SAM of reads aligned to a template, writes a FASTA of reads translated
    within the open reading frame of their alignment.

    The translations are validated based on: 1) Mapping quality, 2) Matching to template FRW4 (within a margin of error), and 3) In-alignment truncations.

    FASTA entries are returned with an ID corresponding to the original read ID. Additionally, the description includes three fields:

    * ``framework``: The framework assigned during alignment by ``minimap2`` (The name given in the template fasta file).
    * ``5-umi``: The 5'-UMI of the read, if found. Uses the pattern ``TTTVVVVTTVVVVTTVVVVTTVVVVTTT``, where V ∈ [A,C,G]
    * ``3-umi``: The 3'-UMI of the read, if found. Uses the pattern ``AAABBBBAABBBBAABBBBAABBBBAAA``, where B ∈ [T,C,G]

    :param input_sam: Path for the input SAM.
    :param output_fasta: Path for the output FASTA. If `None`, reuses basename of ``input_fasta``  (i.e. ``sample1.sam`` => ``sample1.fasta``).
    :param templates_fasta: Path for the templates used in the input SAM alignment.
    :param fasta_type: Whether the templates are written as
    nucleotides or amino acids.
    :param val_min_mapq: The minimum ``minimap2`` mapping quality to accept (Maximum ``60``).
    :param val_n_last: How many of the C-ter template residues to use for FRW4 validation.
    :param val_n_sub: How many substitutions to tolerate for FRW4 validation.
    """

    # Reusing input file name, if no --output-fasta is passed
    if output_fasta is None:
        output_fasta = os.path.basename(input_sam).replace('.sam', '.fasta')

    # Loading patterns for FRW4 validation
    validation_patterns = load_validation_end_patterns(templates_fasta,
                                                       val_n_last, val_n_sub,
                                                       fasta_type)

    # Initializing a FASTA record container
    records = []

    # Looping over aligned reads
    with AlignmentFile(input_sam, 'r') as sam:
        for read in sam:

            # Skipping badly mapped reads
            if not read.is_mapped or read.mapq < val_min_mapq:
                continue

            # TODO: Make separate script for FRW4 validation
            # Skipping reads that do not match template FRW4
            aa_seq = translate_aligned_read(read)
            framework = read.reference_name
            ab_seq = re.search(validation_patterns[framework], aa_seq)
            if ab_seq is None:
                continue
            else:
                ab_seq = Seq(ab_seq.group())

            # Skipping truncated reads
            if '*' in ab_seq:
                continue

            # Formatting translation as FASTA record
            umi_5, umi_3 = umis(str(read.query_sequence),
                                not_found='N/A')

            record = SeqRecord(ab_seq,
                               id=read.query_name,
                               description=f"framework={framework} 5-umi={umi_5} 3-umi={umi_3}")
            records.append(record)

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
            -t/--templates-fasta    TEMPLATES.FASTA
            [-f/--fasta-type        [dna]                    ]
            [-m/--val-min-mapq      [60]                     ]   
            [-n/--val-n-last        [11]                     ]     
            [-s/--val-n-sub         [1]                      ]      
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
    parser.add_argument('-t', '--templates-fasta',
                        type=str,
                        required=True,
                        help="Templates used for read alignment.")
    parser.add_argument('-f', '--fasta-type',
                        type=str,
                        required=False,
                        default='dna',
                        choices=['dna', 'aa'],
                        help="Whether templates are given as nucleotides [dna] or amino acids [aa].")
    parser.add_argument('-m', '--val-min-mapq',
                        type=int,
                        required=False,
                        default=60,
                        help="The minimum minimap2 mapping quality to accept (Maximum 60).")
    parser.add_argument('-n', '--val-n-last',
                        type=int,
                        required=False,
                        default=11,
                        help="How many of the C-ter template residues to use for FRW4 validation.")
    parser.add_argument('-s', '--val-n-sub',
                        type=int,
                        required=False,
                        default=1,
                        help="How many substitutions to tolerate for FRW4 validation.")

    args = parser.parse_args()

    # ············································································· #

    # Translating SAMs sequentially
    translate_sam(**vars(args))
