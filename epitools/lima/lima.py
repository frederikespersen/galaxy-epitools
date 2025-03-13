#!/usr/bin/env python3

# ················································································· #

# Imports
import os
import glob
import shutil
import argparse

# ················································································· #

def process_fastq_files(input_dir: str,
                        fastq_prefix: str,
                        output_dir: str) -> None:
    """
    Processes FASTQ files in the input directory, renames them,
    and copies them to the output directory.

    :param input_dir: Path to the directory containing Lima-demultiplexed FASTQ files.
    :param fastq_prefix: Prefix given by Lima to FASTQ files.
    :param output_dir: Path to the directory to save renamed FASTQ files.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through all FASTQ files in the input directory
    for fastq in glob.glob(os.path.join(input_dir, "*.fastq")):

        # Extract the new file name
        fastq_out = os.path.basename(fastq).replace(fastq_prefix+'.', '').split('--')[0].replace('.fastq','')

        # Construct the output file path
        output_file = os.path.join(output_dir, f"{fastq_out}.fastq")

        # Copy and rename the file
        shutil.copy(fastq, output_file)


# ················································································· #

if __name__ == "__main__":

    # Parsing input arguments
    parser = argparse.ArgumentParser(
        prog="process-fastq",
        usage="""
            python process-fastq.py
                [-h]
                -i/--input-dir         INPUT_DIR
                -o/--output-dir        OUTPUT_DIR
                -p/--fastq-prefix      FASTQ_PREFIX
                """,
        description="Processes and renames FASTQ files based on specific naming rules.",
        epilog="// Frederik Espersen Knudsen, 2024",
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        type=str,
        required=True,
        help="Path to the directory containing Lima-demultiplexed FASTQ files.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Path to the directory to save renamed FASTQ files.",
    )
    parser.add_argument(
        "-p",
        "--fastq-prefix",
        type=str,
        required=True,
        help="Prefix given by Lima to FASTQ files.",
    )

    args = parser.parse_args()

    # ············································································· #

    # Process FASTQ files
    process_fastq_files(input_dir=args.input_dir, output_dir=args.output_dir, fastq_prefix=args.fastq_prefix)
