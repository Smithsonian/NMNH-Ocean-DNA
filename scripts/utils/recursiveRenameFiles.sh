#!/bin/bash
# script to replace string in file names recursively within a specified directory
# author: Dan MacGuigan
# consider making a backup of the directory before running this script
# script takes 2 arguments:
# 1) txt file with two columns, first column is old string in file names, second column is new string in file names
# 2) directory to recursively search

###########
# Help
###########
Help()
{
# display help
echo ''
echo 'script to replace string in file names recursively within a specified directory'
echo 'author: Dan MacGuigan'
echo ''
echo 'PLEASE consider making a backup of the directory before running this script'
echo 'script takes 2 arguments:'
echo '1) txt file with two columns, first column is old string in file names, second column is new string in file names'
echo '2) directory to recursively search'
echo ''
echo 'usage: recursiveRenameFiles.sh NAMES_FILE.txt /PATH/TO/DIRECTORY/TO/SEARCH'
echo ''
}

# Get the options
while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done



shopt -s globstar # allow recursive wildcard **
while IFS= read -r line; 
do
	t=(${line}) # make array from line in txt file
	echo "old string: ${t[0]}"
	echo "new string: ${t[1]}"
	rename ${t[0]} ${t[1]} ${2}/**
done < ${1}
