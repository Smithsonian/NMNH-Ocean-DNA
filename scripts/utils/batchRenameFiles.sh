#!/bin/bash

# script to rename many files in a specified directory
# can replace all or part of the file names
# author: Dan MacGuigan
# PLEASE consider making a backup of the directory before running this script
# script takes 2 arguments:
# 1) tab-delimited txt file with two columns:
#       first column contains the old file name (or old part of a file name)
#       second column contains the new file name (or new part of a file name)
#       example file:
#          sample1.fasta	sample1_newStuff.fasta
#          s1	newStuff 
#       this will rename the file "sample1.fasta" to "sample1_newStuff.fasta"
#       and it will replace "s1" with "newStuff" in all file names
# 2) directory to recursively search
# usage: batchRenameFiles.sh NAMES_FILE.txt /PATH/TO/DIRECTORY/TO/SEARCH

###########
# Help
###########
Help()
{
# display help
echo ''
echo 'Script to rename many files in a specified directory'
echo 'Can replace all or part of the file names'
echo 'Author: Dan MacGuigan'
echo ''
echo 'PLEASE consider making a backup of the target directory before running this script'
echo ""
echo 'This script has two modes:'
echo '"-m full" will change full file names'
echo '"-m partial" will change part of file names'
echo ''
echo 'Syntax: batchRenameFiles.sh [-h|d|t|m]'
echo 'Options:'
echo "h     Print this Help"
echo "d     Directory containing files to rename"
echo "t     Tab-delimited txt file with two columns:"
echo "        First column: old file names or old strings in file names"
echo "        Second column: new file names or new strings in file names"
echo "m     File rename mode [full|partial]"
echo ''
}

# Get the options
while getopts ":hd:t:m:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      d)
        DIR=${OPTARG};;
      t)
        TABLE=${OPTARG};;
      m)
        MODE=${OPTARG};;
      \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

# check variables were set
if [[ $# -eq 0 ]] ; then
    Help
    exit 0
fi
if [ -z ${DIR+x} ]; 
then 
  echo "Error: please specify a directory (-d)"
  exit
fi
if [ -z ${TABLE+x} ]; 
then 
  echo "Error: please specify a renaming table (-t)"
  exit
fi
if [ -z ${MODE+x} ]; 
then 
  echo "Error: please specify a mode (-m [full|partial])"
  exit
fi
# check variables were set properly
if ! [ -d "${DIR}" ]
then
  echo "Error: ${DIR} does not exist"
  exit
fi
if ! [ -f "${TABLE}" ]
then
  echo "Error: ${TABLE} does not exist"
  exit
fi
if ! [[ "${MODE}" =~ ^(full|partial)$ ]]
then
  echo 'Error: Invalid choice for -m, must be "full" or "partial"'
  exit
fi

if [ ${MODE} == "full" ]
then
	echo "Performing full file name replacement"
	echo ""
	while IFS= read -r line; 
	do
		t=(${line}) # make array from line in txt file
		echo "${t[0]} -> ${t[1]}"
		mv ${DIR}/${t[0]} ${DIR}/${t[1]}
	done < ${TABLE}
elif [ ${MODE} == "partial" ]
then
	echo "Performing partial file name replacement"
	while IFS= read -r line; 
	do
		echo ""
		t=(${line}) # make array from line in txt file
		echo "replacing \"${t[0]}\" with \"${t[1]}\" in file names:"
		find ${DIR}/*${t[0]}* -maxdepth 1 -type f -printf "%f\n" > files.temp.txt
		sed "s/${t[0]}/${t[1]}/g" files.temp.txt > filesNew.temp.txt
		paste files.temp.txt filesNew.temp.txt | sed 's/\t/ -> /' 
		rm files.temp.txt filesNew.temp.txt
		rename ${t[0]} ${t[1]} ${DIR}/*
	done < ${TABLE}
fi

echo ""
echo "done renaming files!"