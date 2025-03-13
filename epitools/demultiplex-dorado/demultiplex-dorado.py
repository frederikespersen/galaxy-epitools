#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import os
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# ················································································· #

def dorado_barcode(read: SeqRecord,
                   rg_z_pos: int=3) -> str | bool:
    """
    Extracts the name of identified barcode assigned by Dorado in the RG:Z: field.
    Requires reads from ``dorado basecaller`` with the ``--kit-name <barcode-kit-name>``
    argument specified.

    :param read: A Bio FASTQ SeqRecord of a Dorado read.
    :param rg_z_pos: The 0-indexed space-delimited position of the RG:Z: field in the read header [Default: ``3``] (See ``dorado_rg_z_pos()``).
    :return: The barcode name if found, else ``False``.
    """
    fields = read.description.split('\t')
    rg_z_tag = fields[rg_z_pos]
    if 'barcode' in rg_z_tag:
        return '_'.join(rg_z_tag.split('_')[-2:])
    else:
        return False


def dorado_rg_z_pos(read: SeqRecord) -> int:
    """
    Determines the 0-indexed space-delimited position of the RG:Z: field in a Dorado read's description.
    Returns a ValueError if the RG:Z: field is not found

    :param read: A Bio FASTQ SeqRecord of a Dorado read.
    :return: The 0-indexed position of the RG:Z: field by space delimitation.
    """
    fields = read.description.split('\t')
    for i, field in enumerate(fields):
        if field[:5] == 'RG:Z:':
            return i
    raise ValueError("RG:Z: field not found in read header!")


def fastq_barcode_split(input_fastq: str,
                        output_dir: str,
                        nonbarcoded: str='unbarcoded') -> None:
    """
    Splits FASTQ file by barcode assigned by Dorado.
    New FASTQ files are named according to barcodes.
    Requires reads from ``dorado basecaller`` with the ``--kit-name <barcode-kit-name>``
    argument specified. Extracts barcode with ``dorado_barcode()``.

    :param input_fastq: Path for FASTQ file with reads to split by barcode.
    :param output_dir: The directory in which to save split FASTQ files.
    :param nonbarcoded: The FASTQ filename for reads without an assigned barcode.
    """

    # Collecting reads according to barcode
    barcodes = {}

    # Parsing FASTQ file(s)
    fastq = SeqIO.parse(input_fastq, format='fastq')
    rg_z_pos = None
    for read in fastq:

        if rg_z_pos is None:
            rg_z_pos = dorado_rg_z_pos(read)

        # Assign barcode to read
        barcode = dorado_barcode(read, rg_z_pos)
        if not barcode:
            barcode = nonbarcoded

        # Initialise barcode container
        if barcode not in barcodes:
            barcodes[barcode] = []

        # Assign read to collection
        barcodes[barcode].append(read)

    # Creating output directory
    os.makedirs(output_dir, exist_ok=True)

    # Writing read collections to separate FASTQ files
    for barcode, reads in barcodes.items():
        with open(os.path.join(output_dir, f'{barcode}.fastq'), 'w') as output_fastq:
            SeqIO.write(reads, output_fastq, format='fastq')

# ················································································· #

if __name__ == '__main__':

    # Getting commandline arguments
    parser = argparse.ArgumentParser(
        prog='demulitplex-dorado',
        usage='''
    python demultiplex-dorado.py
        [-h]
        -i/--input-fastq    INPUT.FASTQ
        -o/--output-dir     OUTPUT_DIR
        ''',
        description="Demultiplex FASTQ reads by Dorado-assigned barcodes.",
        epilog='// Frederik Espersen Knudsen, 2024')
    parser.add_argument('-i', '--input-fastq',
                        type=str,
                        required=True,
                        help="FASTQ file with Dorado-assigned barcodes.")
    parser.add_argument('-o', '--output-dir',
                        required=True,
                        type=str,
                        default='demultiplexed',
                        help="Directory to save demultiplexed FASTQ files in.")

    args = parser.parse_args()

    # ············································································· #

    # Splitting FASTQ input by barcodes
    fastq_barcode_split(args.input_fastq, args.output_dir)
