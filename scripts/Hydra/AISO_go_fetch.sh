#!/bin/bash

# ABOUT ########################################################
# author: Dan macguigan
# contact: macguigand@si.edu
#
# Script to generate custom GetOrganlle reference databases
# From https://github.com/o-william-white/go_fetch
#
# Requires a conda environment "go_fetch", which can be created by:
### conda create -n go_fetch -c bioconda getorganelle biopython trf
#
# Also requires the go_fetch.py script to be in your PATH
# To add the script to your bin directory, run the following:
### cd ~
### git clone https://github.com/o-william-white/go_fetch
### mkdir -p ~/bin
### cp -rf ~/go_fetch/go_fetch.ph ~/bin
#
# Example usage:
# [manually modify variables in the INPUTS section]
# run "bash AISO_go_fetch.sh"

# INPUTS ########################################################

CLADE_NAME="cnidaria" # Name of your clade, only used for naming output folder
NCBI_TAXID="6073" # NCBI taxonomy ID of rank to search for sequences
TARGET="mitochondrion" # Target sequence type. Options = [chloroplast, mitochondrion, ribosomal, ribosomal_complete]
MIN=5 # Minimum number of target sequences to download
MAX=50 # Maximum number of target sequences to download

# SCRIPT ########################################################

conda activate go_fetch

DATE=$(date '+%Y-%m-%d')

mkdir ${CLADE_NAME}_${TARGET}_max${MAX}_${DATE}_refDBs
cd ${CLADE_NAME}_${TARGET}_max${MAX}_${DATE}_refDBs

# refseq
go_fetch.py \
   --taxonomy ${NCBI_TAXID}  \
   --target ${TARGET} \
   --db refseq \
   --min ${MIN} --max ${MAX} \
   --output refseq \
   --overwrite \
   --getorganelle \
   --email user_email@example.com

# genbank
go_fetch.py \
   --taxonomy ${NCBI_TAXID}  \
   --target ${TARGET} \
   --db genbank \
   --min ${MIN} --max ${MAX} \
   --output genbank \
   --overwrite \
   --getorganelle \
   --email user_email@example.com

# combine genbank and refseq results
cat refseq/seed.fasta genbank/seed.fasta > seed.fasta
cat refseq/gene.fasta genbank/gene.fasta > gene.fasta

cd ..

# OUTPUT ########################################################
echo ""
echo "Results saved to ./${CLADE_NAME}_${TARGET}_max${MAX}_${DATE}_refDBs"
echo ""
echo "If you are using AISO_GetOrganelle.sh, include the following lines in your INPUTS section of the script:"
echo ""
echo "SEED_DB=\"${PWD}/${CLADE_NAME}_${TARGET}_max${MAX}_${DATE}_refDBs/seed.fasta\""
echo "LABEL_DB=\"${PWD}/${CLADE_NAME}_${TARGET}_max${MAX}_${DATE}_refDBs/gene.fasta\""
echo ""
