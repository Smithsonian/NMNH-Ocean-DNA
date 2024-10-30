#!/bin/bash
#$ -N GetOrganelle
#$ -o system.out
#$ -e error.out
#$ -j y
#$ -terse
#$ -notify
#$ -q sThC.q
#$ -pe mthread 6
#$ -l mres=12G,h_data=2G,h_vmem=2G 
#$ -S /bin/bash
#$ -cwd

## ABOUT ################################################
# author: Dan macguigan
# contact: macguigand@si.edu
#
# Script to run GetOrganelle on NMNH Hydra cluster,
# adapted from MitoPilot. Assumes you have GetOrganlle 
# installed and in your PATH (i.e. you can call 
# get_organelle_from_reads.py from anywhere)
# 
# To install GetOrganelle on Hydra:
# conda create --name GetOrganelle bioconda::getorganelle
# conda activate GetOrganelle
# get_organelle_from_reads.py -h 
# 
# More GetOrganelle details here:
# https://github.com/Kinggerm/GetOrganelle
#
# Example usage:
# [manually modify variables in the INPUTS section] 
# qsub AISO_GetOrganelle_Hydra.sh
#
## INPUTS ################################################

# main working directory, output will be written here
WORD_DIR="/home/macguigand/pool/GetOrganelle_test"

# sample info
SAMPLE_ID="SRR22396940" # could be anything, no spaces
R1="/pool/public/genomics/macguigand/GetOrganelle_test/SRR22396940_preprocess_R1.fastq.gz"
R2="/pool/public/genomics/macguigand/GetOrganelle_test/SRR22396940_preprocess_R2.fastq.gz"

# reference databases, see link below for more details
# https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference
SEED_DB="/full/path/to/your/seed_database.fasta"
LABEL_DB="/full/path/to/your/genes_database.fasta"

THREADS=6 # number of threads/CPUs, must match "-pe mthread" in the UGE block above


## SCRIPT ################################################

cd ${WD}

conda activate GetOrganelle # assumes you created a conda environment called "GetOrganelle" for install

mkdir -p "${SAMPLE_ID}/assemble"

# run GetOrganelle, default settings from MitoPilot
get_organelle_from_reads.py \
 -1 ${R1} \
 -2 ${R2} \
 -o ${SAMPLE_ID}/assemble/ \
 --overwrite \
 -s ${SEED_DB} \
 --genes ${LABEL_DB} \
 -t ${THREADS} \
 -F 'anonym' \
 -R 10 \
 -k '21,45,65,85,105,115' \
 --larger-auto-ws \
 --expected-max-size 20000 \
 --target-genome-size 16500
 
# summarize results
summary_get_organelle_output.py ${SAMPLE_ID}/assemble -o ${SAMPLE_ID}/assemble/${SAMPLE_ID}_summary.txt


## OUTPUT ################################################
echo ""
echo "GetOrganlle COMPLETE!"
echo "assembly located in ${SAMPLE_ID}/assemble"
echo "see ${SAMPLE_ID}/assemble/${SAMPLE_ID}_summary.txt for stats"
