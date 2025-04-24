#!/bin/bash
#$ -N GetOrganelle
#$ -o GetOrganelle.out
#$ -e GetOrganelle.err
#$ -j y
#$ -terse
#$ -notify
#$ -pe mthread 6
#$ -q sThM.q
#$ -l mres=432G,h_data=72G,h_vmem=72G,himem
#$ -S /bin/bash
#$ -cwd

# ABOUT ################################################
# author: Dan MacGuigan
# contact: macguigand@si.edu
#
# Script to run GetOrganelle on NMNH Hydra cluster,
# adapted from MitoPilot. Assumes you have GetOrganlle 
# installed and in your PATH (i.e. you can call 
# get_organelle_from_reads.py from anywhere)
# 
# To install GetOrganelle on Hydra:
### conda create --name GetOrganelle bioconda::getorganelle
### conda activate GetOrganelle
### get_organelle_from_reads.py -h 
# 
# More GetOrganelle details here:
# https://github.com/Kinggerm/GetOrganelle
#
# Example usage:
# [manually modify variables in the INPUTS section] 
# run "qsub AISO_GetOrganelle_Hydra.sh"

# INPUTS ################################################

# main working directory, output will be written here
WD="/pool/public/genomics/macguigand/GetOrganelle_invert_test/Cnidaria/Cnidaria-Pandea-hydroids-1525150"

# sample info
SAMPLE_ID="Cnidaria-Pandea-hydroids-1525150" # could be anything, no spaces
R1="/pool/public/genomics/macguigand/GetOrganelle_invert_test/Cnidaria/Cnidaria-Pandea-hydroids-1525150/rawData/Cnidaria-Pandea-hydroids-1525150_S20_R1_001.fastq.gz"
R2="/pool/public/genomics/macguigand/GetOrganelle_invert_test/Cnidaria/Cnidaria-Pandea-hydroids-1525150/rawData/Cnidaria-Pandea-hydroids-1525150_S20_R2_001.fastq.gz"

# reference databases, see link below for more details
# https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference
SEED_DB="/home/macguigand/pool/GetOrganelle_invert_test/Cnidaria/Cnidaria-Pandea-hydroids-1525150/cnidaria_mitochondrion_max50_2024-11-04_refDBs/seed.fasta"
LABEL_DB="/home/macguigand/pool/GetOrganelle_invert_test/Cnidaria/Cnidaria-Pandea-hydroids-1525150/cnidaria_mitochondrion_max50_2024-11-04_refDBs/gene.fasta"

THREADS=6 # number of threads/CPUs, must match "-pe mthread" in the UGE block above


# SCRIPT ################################################

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
 --target-genome-size 16500 \
 --which-spades /home/macguigand/SPAdes-4.0.0-Linux/bin # NEED TO FIX
 
# summarize results
summary_get_organelle_output.py ${SAMPLE_ID}/assemble -o ${SAMPLE_ID}/assemble/${SAMPLE_ID}_summary.txt


# OUTPUT ################################################
echo ""
echo "GetOrganlle COMPLETE!"
echo "assembly located in ${SAMPLE_ID}/assemble"
echo "see ${SAMPLE_ID}/assemble/${SAMPLE_ID}_summary.txt for stats"
