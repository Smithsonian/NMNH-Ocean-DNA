#!/bin/bash
#$ -N BLAST_18s-28s
#$ -o BLAST_18s-28s_loop.out
#$ -e BLAST_18s-28s_loop.err
#$ -j y
#$ -terse
#$ -notify
#$ -q sThM.q
#$ -l mres=32G,h_data=32G,h_vmem=32G,himem
#$ -S /bin/bash
#$ -cwd

# ABOUT ################################################
# script to loop through all contigs/scaffolds in a directory
# run BLAST for 18s and 28s sequences
# and extract and rename hits
#
# require the BLAST_18s-28s_hydra.sh
# <<LINK>>
#
# author: Dan MacGuigan
# contact: macguigand@si.edu

# INPUTS ################################################

# working directory containing contig/scaffold FASTA files
DIR="/pool/public/genomics/macguigand/BLAST_testing/scaffolds"

# FASTA file suffix (e.g. "fasta", "fa", "fas")
# must be the same for all files in DIR
SUFFIX="fasta"

# full path to 18s and 28s query sequences
rRNA_S="/pool/public/genomics/macguigand/BLAST_testing/18s.fasta"
rRNA_L="/pool/public/genomics/macguigand/BLAST_testing/28s.fasta"

# full path to BLAST_18s-28s_hydra.sh script
BLAST_SCRIPT="/pool/public/genomics/macguigand/BLAST_testing/BLAST_18s-28s_hydra.sh"

# SCRIPT ################################################
# DO NOT MODIFY #########################################

# loop through every file ending with SUFFIX in DIR
cd ${DIR}
for FILE in *${SUFFIX}
do
    ID=${FILE%".$SUFFIX"} # sample ID minus FASTA suffix
    echo ""
    echo "bash ${BLAST_SCRIPT} -c ${FILE} -i ${ID} -s ${rRNA_S} -l ${rRNA_L}"
    bash ${BLAST_SCRIPT} -c ${FILE} -i ${ID} -s ${rRNA_S} -l ${rRNA_L}
done

# OUTPUT ################################################
# DO NOT MODIFY #########################################

echo ""
echo "DONE processing all samples"
echo "BLAST hits written to:"
echo "${DIR}/BLAST_hits"