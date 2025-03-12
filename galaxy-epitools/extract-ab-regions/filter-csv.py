#!usr/bin/env python3

# ················································································· #

# Imports
import argparse
import pandas as pd

# ················································································· #

def filter_csv(input_file: str, filters: dict) -> pd.DataFrame:
    """Filters a CSV file by column values."""
    df = pd.read_csv(input_file)
    
    for column, value in filters.items():
        df = df[df[column] == value]
    
    return df

# ················································································· #

if __name__ == '__main__':

    # Parsing input arguments
    parser = argparse.ArgumentParser(
        prog='filter-csv',
        usage='''
            python filter_csv.py
                [-h]
                -i/--input        INPUT.[CSV]
                COLUMN=VALUE ...
            ''',
        description="Filter a CSV file by column values.",
        epilog='// Frederik Espersen Knudsen, 2025'
    )

    parser.add_argument('-i', '--input',
                        type=str,
                        required=True,
                        help="Path to the input CSV file.")
    parser.add_argument('filters',
                        nargs='+',
                        help="Column filters in the format 'Column=Value'.")

    args = parser.parse_args()
    
    filters = {}
    for f in args.filters:
        column, value = f.split("=", 1)
        filters[column] = value
    
    # Filtering and printing output
    filtered_df = filter_csv(args.input, filters)
    filtered_df.to_csv(args.input, index=False)