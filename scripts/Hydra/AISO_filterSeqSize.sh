#!/bin/bash

# ABOUT ################################################
# author: Dan MacGuigan
# contact: macguigand@si.edu
#
# Script to filter FASTA sequences by minimum size.
# Will run on all FASTA files in the supplied directory
# and write results to new files.
# 
# Written specifically for NMNH Hydra cluster, but can
# work on any computer with seqtk installed.
# https://github.com/lh3/seqtk
#

# INPUTS ################################################

# directory containing FASTA files you wish to filter by sequence size
DIR="/pool/public/genomics/macguigand/test/files"

# file suffix for your FASTA sequence files in DIR, must include the period
SUFFIX=".fasta"

# minumum sequence legth (in bp)
MIN_L=1000

# SCRIPT ################################################

module load bioinformatics/seqtk/1.4

mkdir ${DIR}/L${MIN_L}

echo "Filtering files:"

for file in "${DIR}"/*"${SUFFIX}"
do
  echo ${file}
  base_name=$(basename ${file})
  file_prefix=${base_name%"$SUFFIX"} # get base file name
  seqtk seq -L "${MIN_L}" ${file} > ${DIR}/L${MIN_L}/${file_prefix}_L${MIN_L}${SUFFIX} # filter
done

# OUTPUT ################################################
echo ""
echo "Sequencing filtering complete!"
echo "Results written to"
echo "${DIR}/L${MIN_L}"
echo "with file endings \"_L${MIN_L}${SUFFIX}\""

