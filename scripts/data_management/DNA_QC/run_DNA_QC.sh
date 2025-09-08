#!/bin/bash
##################################################################
## ABOUT
# script to generate QC report for FASTQ DNA sequence data
# requires Nextflow to be installed
# see instructions here: 

# NOTE: if you have many samples, consider running this as a batch job
# you will need to use the special workflow manager queue "lTWFM.sq"

##################################################################
## INPUTS
# directory containing paired DNA sequence reads
# files must be FASTQ and gzipped (i.e. end in ".fastq.gz")
RAW_DATA_DIR="/pool/public/genomics/macguigand/QC_test/raw_data"

# suffix for your paired forward and reverse reads
# if your paired reads are "sampleA_R1.fastq.gz" and "sampleA_R2.fastq.gz" 
# the suffix below will work
READS_SUFFIX="_R{1,2}.fastq.gz"

##################################################################
# SCRIPT
source ~/.bashrc
module load tools/java/21.0.2 # required for Nextflow
nextflow run ~/DNA_QC_script/DNA_QC.nf -entry QC -c ~/DNA_QC_script/Hydra.nf.config --raw_data_dir ${RAW_DATA_DIR} --reads_suffix ${READS_SUFFIX}
