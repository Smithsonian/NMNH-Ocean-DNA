#!/bin/bash

# ABOUT ################################################
# author: Dan MacGuigan
# contact: macguigand@si.edu
# Script to find and rename contigs/scaffolds containing 18s and/or 28s genes
# Run "bash BLAST_18s-28s_hydra.sh -h" for help

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo ""
   echo "Script to find and rename contigs/scaffolds containing"
   echo "18s and/or 28s genes"
   echo ""
   echo "author: Dan MacGuigan"
   echo "contact: macguigand@si.edu"
   echo ""
   echo "Options:"
   echo "c   FASTA file containing genomics scaffolds or contigs"
   echo "i   sample ID, will be used to name resulting files and BLAST hits"
   echo "s   FASTA file containing query 18s sequence"
   echo "l   FASTA file containing query 28s sequence"
   echo "h   Print this Help"
   echo ""
   echo "Usage:"
   echo "bash BLAST_18s-28s_hydra.sh -c my_contigs.fasta -i my_sample_ID -s my_18s.fasta -l my_28s.fasta"
}

############################################################
# Process the input options.                               #
############################################################
# Get the options

unset -v CONTIGS
unset -v SMALL
unset -v LARGE

while getopts "hc:s:i:l:c:" option; do
   case ${option} in
      h) # display Help
         Help
         exit;;
      c) # contig FASTA
         CONTIGS=$OPTARG
         ;;
      i) # sample ID
         ID=$OPTARG
         ;;
      s) # small subunit (18s) FASTA
         SMALL=$OPTARG
         ;;
      l) # large subunit (28s) FASTA
         LARGE=$OPTARG
         ;;
      \?) # Invalid option
        echo "Error: Invalid option"
        Help
        exit;;
   esac
done

shift "$(( OPTIND - 1 ))"

if [ -z "$CONTIGS" ] || [ -z "$SMALL" ] || [ -z "$LARGE" ] || [ -z "$ID" ]; then
        echo 'WARNING: Missing required arguments' >&2
        Help
        exit 1
fi

# SCRIPT ################################################

echo "PROCESSING sample: ${CONTIGS}"

# load BLAST module on Hydra
module load bio/blast/2.15.0

# create BLAST database from genome assembly
makeblastdb -in "${CONTIGS}" -dbtype "nucl" -out "temp_blast_db"

# BLAST 18s
mkdir -p 18s_BLAST_results
blastn -query "${SMALL}" -db "temp_blast_db" -outfmt 6 -out "18s_BLAST_results/${CONTIGS}.18s.blast.txt"

# BLAST 28s
mkdir -p 28s_BLAST_results
blastn -query "${LARGE}" -db "temp_blast_db" -outfmt 6 -out "28s_BLAST_results/${CONTIGS}.28s.blast.txt"

# list of BLAST hits
mkdir -p blast_temp
cut "18s_BLAST_results/${CONTIGS}.18s.blast.txt" -f2 | sort | uniq > blast_temp/18s_temp.txt
cut "28s_BLAST_results/${CONTIGS}.28s.blast.txt" -f2 | sort | uniq > blast_temp/28s_temp.txt

# join lists
comm -12 blast_temp/18s_temp.txt blast_temp/28s_temp.txt  > blast_temp/18s-28s_contigs.txt
comm -23 blast_temp/18s_temp.txt blast_temp/28s_temp.txt  > blast_temp/18s_contigs.txt
comm -13 blast_temp/18s_temp.txt blast_temp/28s_temp.txt  > blast_temp/28s_contigs.txt

# extract sequences from CONTIGS
module load bio/seqtk/1.4
seqtk subseq "${CONTIGS}" blast_temp/18s_contigs.txt > blast_temp/18s_temp.fasta
seqtk subseq "${CONTIGS}" blast_temp/28s_contigs.txt > blast_temp/28s_temp.fasta
seqtk subseq "${CONTIGS}" blast_temp/18s-28s_contigs.txt > blast_temp/18s-28s_temp.fasta

# modify contig names
sed -i "s/>/>$ID 18s_blast_hit /g" blast_temp/18s_temp.fasta
sed -i "s/>/>$ID 28s_blast_hit /g" blast_temp/28s_temp.fasta
sed -i "s/>/>$ID 18s-28s_blast_hit /g" blast_temp/18s-28s_temp.fasta

# merge extracted seqs and write to output directory
mkdir -p BLAST_hits
cat blast_temp/18s-28s_temp.fasta blast_temp/18s_temp.fasta blast_temp/28s_temp.fasta > BLAST_hits/${ID}.blastHits.fasta
 
# cleanup
rm -rf blast_temp
rm temp_blast_db.*
rm -rf *_BLAST_results

# OUTPUTS ################################################

echo ""
echo "FINISHED sample: ${CONTIGS}"
echo "results written to:"
echo "${PWD}/BLAST_hits/${ID}.blastHits.fasta"