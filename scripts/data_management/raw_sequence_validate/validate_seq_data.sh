#!/bin/bash

# ABOUT ################################################
# author: Dan MacGuigan
# contact: macguigand@si.edu
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
# bash validate_seq_data.sh map_file.txt

# HELP FUNCTION ########################################
# This function displays the help message for the script.
show_help() {
cat << EOF
Usage: ${0##*/} [-h] MAP_FILE

DESCRIPTION:
    Validates Ocean DNA raw sequence data against its corresponding metadata.
    This script checks for correct naming, file existence, and ensures that
    the contents of the metadata file and sequence data directory match.

ARGUMENTS:
  MAP_FILE
      Path to a text file mapping metadata to sequence data directories.
      It requires a header row (which is ignored). Each subsequent line
      should contain two space-separated columns:

      Column 1: The metadata CSV filename.
                - Must be in "/store/public/oceandna/raw_sequence_metadata"
                - Must end with "_mapfile.csv"
                - Header must start with "ID,R1,R2,Taxon,UniqueID"

      Column 2: The raw sequence data directory name.
                - Must be in "/store/public/oceandna/raw_sequence_data"
                - Must contain one or more .fastq.gz files

OPTIONS:
  -h, --help
      Display this help message and exit.

EXAMPLE:
    bash ${0##*/} my_projects.txt
EOF
}

# INPUTS ################################################

# Parse command-line options for help flag
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# Check if a map file argument was provided
if [ -z "$1" ]; then
    echo "Error: No map file provided." >&2
    show_help
    exit 1
fi

MAP_FILE=$1

# Check if the map file exists and is readable
if [ ! -f "$MAP_FILE" ]; then
    echo "Error: Map file '$MAP_FILE' not found." >&2
    exit 1
fi

# SCRIPT ################################################

N1=0
while read LINE;
do
    # Skip empty lines
    if [ -z "$LINE" ]; then
        continue
    fi

    if (( N1 > 0 )); then # skip first line containing header
        WARN_COUNT=0
        DATA=($LINE)
        METADATA=${DATA[0]}
        SEQDATA=${DATA[1]}

        echo ""
        echo "VALIDATING"
        echo "metadata: /store/public/oceandna/raw_sequence_metadata/${METADATA}" # metadata
        echo "raw data directory: /store/public/oceandna/raw_sequence_data/${SEQDATA}" # raw data directory

        # cleanup old reports
        PRJ=${METADATA%"_mapfile.csv"}
        rm -f "${PRJ}.report.txt"

        # check metadata CSV suffix
        if [[ "$METADATA" == *"_mapfile.csv" ]]; then
            echo "✅ $METADATA ends in _mapfile.csv"
        else
            echo "❌ $METADATA does not end in _mapfile.csv"
            echo "WARNING: $METADATA does not end in _mapfile.csv" >> ${PRJ}.report.txt
            continue
        fi

        #### SEQUENCE DATA VALIDATION
        # check if sequence directory exists
        if [ -d /store/public/oceandna/raw_sequence_data/${SEQDATA} ]; then
            echo "✅ /store/public/oceandna/raw_sequence_data/${SEQDATA} exists"
        else 
            echo "❌ /store/public/oceandna/raw_sequence_data/${SEQDATA} does not exist"
            echo "WARNING: /store/public/oceandna/raw_sequence_data/${SEQDATA} does not exist" >> ${PRJ}.report.txt
            continue
        fi
        # check if sequence directory contains .fastq.gz files

        fastq_files=$(find "/store/public/oceandna/raw_sequence_data/${SEQDATA}" -name "*.fastq.gz" -type f 2>/dev/null)
        if [ -n "$fastq_files" ]; then
            echo "✅ /store/public/oceandna/raw_sequence_data/${SEQDATA} contains .fastq.gz files"
        else
            echo "❌ /store/public/oceandna/raw_sequence_data/${SEQDATA} contains no .fastq.gz files"
            echo "WARNING: /store/public/oceandna/raw_sequence_data/${SEQDATA} contains no .fastq.gz files" >> ${PRJ}.report.txt
            continue
        fi

        # loop through each .fastq.gz file in raw data directory and see if it is present in the metadata CSV
        for FILE in /store/public/oceandna/raw_sequence_data/${SEQDATA}/*.fastq.gz
        do
            FILE_BASE=$(basename $FILE)
            if grep -q ${FILE_BASE} /store/public/oceandna/raw_sequence_metadata/${METADATA}; then
                continue
            else
                echo "WARNING: ${FILE_BASE} not in metadata CSV" >> ${PRJ}.report.txt
                WARN_COUNT=$((WARN_COUNT+1))
            fi
        done

        #### METADATA SPREADSHEET VALIDATION

        # check that first five columns of metadata CSV are "ID", "R1", "R2", "Taxon", "UniqueID"
        header_line=$(head -n 1 "/store/public/oceandna/raw_sequence_metadata/${METADATA}")
        if [ -z "$header_line" ]; then
            echo "❌ Metadata file appears to be empty"
            echo "ERROR: Metadata file appears to be empty" >> ${PRJ}.report.txt
            continue
        fi
        if [[ "$header_line" == "ID,R1,R2,Taxon,UniqueID"* ]]; then
            echo "✅ first five columns of metadata CSV are ID, R1, R2, Taxon, and UniqueID"
        else
            echo "❌ first five columns of metadata CSV are not ID, R1, R2, Taxon, and UniqueID"
            echo "ERROR: first five columns of metadata CSV are not ID, R1, R2, Taxon, and UniqueID" >> ${PRJ}.report.txt
            continue
        fi 

        # loop through metadata CSV and check if R1 and R2 exist in the raw data directory
        N2=0
        while IFS= read -r L; # Use IFS= and -r for safer reading
        do
            if (( N2 > 0 )); then # skip first line containing header
                # Ensure line is not empty before processing
                if [ -n "$L" ]; then
                    IFS=',' read -r -a SAMPLE <<< "$L"
                    if [ ! -f /store/public/oceandna/raw_sequence_data/${SEQDATA}/${SAMPLE[1]} ]; then
                        echo "WARNING: ${SAMPLE[1]} not found in raw data directory"  >> ${PRJ}.report.txt
                        WARN_COUNT=$((WARN_COUNT+1))
                    fi
                    if [ ! -f /store/public/oceandna/raw_sequence_data/${SEQDATA}/${SAMPLE[2]} ]; then
                        echo "WARNING: ${SAMPLE[2]} not found in raw data directory"  >> ${PRJ}.report.txt
                        WARN_COUNT=$((WARN_COUNT+1))
                    fi
                fi
            fi
            N2=$((N2+1))
        done < /store/public/oceandna/raw_sequence_metadata/${METADATA}
        
        if (( WARN_COUNT > 0 )); then
            echo "⚠️  ${WARN_COUNT} problem(s) found"
            echo "see ${PRJ}.report.txt for details"
        else
            echo "✅ no problems to report"
            echo "no problems to report" >> ${PRJ}.report.txt
        fi
        echo "COMPLETE!"
    
    fi
    N1=$((N1+1))
done < ${MAP_FILE}

# OUTPUT ################################################

echo ""
echo "FINISHED validating all projects in ${MAP_FILE}"