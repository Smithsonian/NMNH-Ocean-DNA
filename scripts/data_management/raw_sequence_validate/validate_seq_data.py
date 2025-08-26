#!/usr/bin/env python3

# ABOUT ################################################
# author: Dan MacGuigan (Bash script), converted to Python by Gemini Coding Partner
#
# script to validate Ocean DNA raw sequence data
#
# requires a TXT file with two columns
#   - metadata: CSV metadata file within "/store/public/oceandna/raw_sequence_metadata".
#               CSV must have header row with first five columns "ID", "R1", "R2", "Taxon", and "UniqueID" (case sensitive).
#               File must end in "_mapfile.csv" (case sensitive).
#   - datadir: directory within "/store/public/oceandna/raw_sequence_data" containing .fastq.gz files
#
# Example usage:
# python validate_seq_data.py map_file.txt
########################################################

import sys
import os
import glob
import csv
import argparse

# Define the base paths for data and metadata.
# Modify these if your directory structure is different.
METADATA_BASE_DIR = "/store/public/oceandna/raw_sequence_metadata"
SEQDATA_BASE_DIR = "/store/public/oceandna/raw_sequence_data"

def write_report(file_path, messages):
    """Writes a list of messages to the specified report file."""
    try:
        with open(file_path, 'w') as f:
            for msg in messages:
                f.write(msg + "\n")
    except IOError as e:
        print(f"Error: Could not write to report file {file_path}. Reason: {e}")

def validate_project(metadata_filename, seq_data_dir_name):
    """
    Runs all validation checks for a single project entry.
    - metadata_filename (str): The name of the metadata CSV file.
    - seq_data_dir_name (str): The name of the sequence data directory.
    """
    print("\nVALIDATING")
    metadata_path = os.path.join(METADATA_BASE_DIR, metadata_filename)
    seq_data_path = os.path.join(SEQDATA_BASE_DIR, seq_data_dir_name)

    print(f"metadata: {metadata_path}")
    print(f"raw data directory: {seq_data_path}")

    # Derive the project name and the path for the report file.
    prj_name = metadata_filename.removesuffix("_mapfile.csv")
    report_file_path = f"{prj_name}.report.txt"

    # A list to hold all error and warning messages for the report.
    report_messages = []

    # --- CRITICAL CHECKS ---
    # If any of these fail, we stop validating this project and move to the next.

    # 1. Check metadata CSV suffix
    if not metadata_filename.endswith("_mapfile.csv"):
        msg = f"WARNING: {metadata_filename} does not end in _mapfile.csv"
        print(f"❌ {metadata_filename} does not end in _mapfile.csv")
        write_report(report_file_path, [msg])
        return

    print(f"✅ {metadata_filename} ends in _mapfile.csv")

    # 2. Check if sequence directory exists
    if not os.path.isdir(seq_data_path):
        msg = f"WARNING: {seq_data_path} does not exist"
        print(f"❌ {seq_data_path} does not exist")
        write_report(report_file_path, [msg])
        return

    print(f"✅ {seq_data_path} exists")

    # 3. Check if sequence directory contains .fastq.gz files
    # The glob module finds all pathnames matching a specified pattern.
    fastq_files_in_dir = glob.glob(os.path.join(seq_data_path, "*.fastq.gz"))
    if not fastq_files_in_dir:
        msg = f"WARNING: {seq_data_path} contains no .fastq.gz files"
        print(f"❌ {seq_data_path} contains no .fastq.gz files")
        write_report(report_file_path, [msg])
        return

    print(f"✅ {seq_data_path} contains .fastq.gz files")

    # 4. Check metadata file and header
    expected_header_start = "ID,R1,R2,Taxon,UniqueID"
    try:
        with open(metadata_path, 'r', newline='') as f:
            header = f.readline().strip()
            if not header:
                msg = "ERROR: Metadata file appears to be empty"
                print("❌ Metadata file appears to be empty")
                write_report(report_file_path, [msg])
                return

            if not header.startswith(expected_header_start):
                msg = f"ERROR: first five columns of metadata CSV are not {expected_header_start}"
                print(f"❌ first five columns of metadata CSV are not {expected_header_start}")
                write_report(report_file_path, [msg])
                return

            print(f"✅ first five columns of metadata CSV are {expected_header_start}")

    except FileNotFoundError:
        msg = f"ERROR: Metadata file not found at {metadata_path}"
        print(f"❌ Metadata file not found at {metadata_path}")
        write_report(report_file_path, [msg])
        return

    # --- NON-CRITICAL CHECKS ---
    # These checks add warnings to the report but don't stop the validation process.

    # 5. Check if every .fastq.gz file is listed in the metadata
    with open(metadata_path, 'r') as f:
        metadata_content = f.read()

    for file_path in fastq_files_in_dir:
        file_base = os.path.basename(file_path)
        if file_base not in metadata_content:
            report_messages.append(f"WARNING: {file_base} not in metadata CSV")

    # 6. Check if R1 and R2 files from metadata exist in the sequence directory
    with open(metadata_path, 'r', newline='') as f:
        reader = csv.reader(f)
        next(reader)  # Skip the header row
        for i, row in enumerate(reader, start=2): # Start from line 2 for logging
            if len(row) < 3:
                report_messages.append(f"WARNING: Malformed row {i} in metadata: {','.join(row)}")
                continue

            # row[1] is the 'R1' column, row[2] is the 'R2' column
            r1_file, r2_file = row[1], row[2]

            if not os.path.isfile(os.path.join(seq_data_path, r1_file)):
                report_messages.append(f"WARNING: {r1_file} (from metadata) not found in raw data directory")

            if not os.path.isfile(os.path.join(seq_data_path, r2_file)):
                report_messages.append(f"WARNING: {r2_file} (from metadata) not found in raw data directory")

    # --- FINAL REPORTING ---
    if report_messages:
        print(f"⚠️  {len(report_messages)} problem(s) found")
        print(f"see {report_file_path} for details")
    else:
        # If no issues were found, add a success message to the report.
        report_messages.append("✅ no problems to report")

    write_report(report_file_path, report_messages)
    print("COMPLETE!")

def main():
    """Main function to parse the map file and trigger validations."""
    # --- NEW: Argument Parsing using argparse ---
    parser = argparse.ArgumentParser(
        description="Validates Ocean DNA raw sequence data against its corresponding metadata.",
        epilog="Example: python validate_seq_data.py my_projects.txt",
        # Keep the formatting of the help text
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Define the positional argument for the map file
    parser.add_argument(
        'map_file',
        metavar='MAP_FILE',
        type=str,
    help="""Path to a text file that maps metadata to sequence data directories.
The file requires a header row (which is ignored). Each following
line must contain two space-separated columns:

  Column 1: The metadata CSV filename.
            - Must be in "/store/public/oceandna/raw_sequence_metadata"
            - Must end with "_mapfile.csv"
            - Header must start with "ID,R1,R2,Taxon,UniqueID"

  Column 2: The raw sequence data directory name.
            - Must be in "/store/public/oceandna/raw_sequence_data"
            - Must contain .fastq.gz files
"""
    )
    
    # Parse the arguments provided by the user
    args = parser.parse_args()
    map_file_path = args.map_file
    # --- END of new argparse section ---

    # Check if the map file exists and is readable
    if not os.path.isfile(map_file_path):
        print(f"Error: Map file not found at {map_file_path}")
        sys.exit(1)

    try:
        with open(map_file_path, 'r') as map_file:
            next(map_file)  # Skip header line of the map file
            for line in map_file:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                parts = line.split()
                if len(parts) < 2:
                    print(f"Skipping malformed line in map file: '{line}'")
                    continue

                metadata_file, seq_data_dir = parts[0], parts[1]
                validate_project(metadata_file, seq_data_dir)

    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
        sys.exit(1)

    print(f"\nFINISHED validating all projects in {map_file_path}")

# This standard Python construct ensures that main() is called
# only when the script is executed directly.
if __name__ == "__main__":
    main()