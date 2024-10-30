#!/bin/bash
# script to replace string in files recursively within a specified directory
# author: Dan MacGuigan
# consider making a backup of the directory before running this script
# script takes 2 arguments:
# 1) txt file with two columns, first column is old string, second column is new string
# 2) directory to recursively search
#    **** TEXT FILE (OPTION 1) MUST NOT BE WITHIN THE DIRECTORY TO SEARCH (OPTION 2)

###########
# Help
###########
Help()
{
# display help
echo ''
echo 'Script to replace strings in files recursively within a specified directory'
echo 'author: Dan MacGuigan'
echo ''
echo '**** I strongly recommend making a backup of the target directory before running this script'
echo ''
echo 'script takes 2 arguments:'
echo '1) txt file with two columns, first column is old string, second column is new string'
echo '2) directory to recursively search'
echo '   **** TEXT FILE (OPTION 1) MUST NOT BE WITHIN THE DIRECTORY TO SEARCH (OPTION 2)'
echo '   **** DO NOT USE THIS SCRIPT TO REPLACE NESTED STRINGS, IT WILL NOT WORK'
echo '   **** for example, you cannot use this script to replace'
echo '   ****    "hell" with "heaven"'
echo '   ****           AND ' 
echo '   ****   "hello" with "goodbye"'
}

# Get the options
while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done

while IFS= read -r line; 
do
	t=(${line}) # make array from line in txt file
	export OLD="${t[0]//./\\.}" # escape period character
	echo "old string: ${OLD}"
	export NEW="${t[1]//./\\.}" # escape period character
	echo "new string: ${NEW}"
        #find ${2} -type f -name "*" -print0 | xargs -0 sed -i -e "s/\b${t[0]}\b/${t[1]}/g"
	perl -p -i -e 's/\b$ENV{OLD}\b/$ENV{NEW}/g' ${2}/*
done < ${1}
