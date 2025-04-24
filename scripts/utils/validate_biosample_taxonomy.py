# ABOUT ########################################################
# author: Dan MacGuigan
# contact: macguigand@si.edu
#
# Script to check taxonomy between all BioSamples and associated
# GenBank nucleotide records in a specific BioProject, then writes
# mismatches to a CSV file.
#
# NOTE: this script can only process a maximum of 10,000 GenBank records
#
# Dependencies: 
#   - python 3 (https://www.python.org/downloads/)
#   - biopython (https://biopython.org/wiki/Download)
#
# This script takes 3 arguments
#   1) a NCBI BioProject number (e.g. PRJNA720393)
#   2) your email address, required for Entrez tools
#   3) your NCBI API key, see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
#
# Example usage:
# python validate_biosample_taxonomy.py "PRJNXXXXX" "YOUR_EMAIL" "YOUR_NCBI_API_KEY"
################################################################

# load libraries
from Bio import SeqIO
from Bio import Entrez
import sys
import fnmatch
from xml.etree import ElementTree as ET
import csv

# command line arguments
bioproject = sys.argv[1]
Entrez.email = sys.argv[2]
Entrez.api_key = sys.argv[3]

# get BioProject stats
count_handle = Entrez.esearch(db="nucleotide", term=f"{bioproject}[BioProject]", retmax=0)
count_results = Entrez.read(count_handle)
count_handle.close()

# Get the count of nucleotide records
n_records = count_results["Count"]

if n_records == 10000:
    print(f"WARNING: BioProject {bioproject} contains {n_records} GenBank records")
    print("          this script can only process a maximum of 10,000 records")
else: 
    print(f"Processing {n_records} GenBank records from BioProject {bioproject}")
    print("this may take a few minutes...")
    print("")

# run esearch to find nucleotides in the BioProject
search_handle = Entrez.esearch(db="nucleotide", term=f"{bioproject}[BioProject]", retmax = 10000)
search_results = Entrez.read(search_handle)
search_handle.close()


# download nucleotide records
id_list = search_results["IdList"]
fetch_handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")

# vars to store results
mismatches = 0
total = 0
no_biosample = 0
nuc_list = []
nuc_tax_list = []
biosample_list = []
biosample_tax_list = []

# parse records
# print biosamples associated with records
for record in SeqIO.parse(fetch_handle, "genbank"):
    try:
        temp = fnmatch.filter(record.dbxrefs, "*BioSample*")[0] # get biosample
    except:
        #print("WARNING - skipping " + record.name + ", no associated BioSample")
        with open('GenBank_records_no_BioSample.txt', 'a') as no_biosample_file:
            no_biosample_file.write(record.name + "\n")
        no_biosample += 1
        continue
    biosample = temp.replace("BioSample:", "") # clean up so only record number
    # download and read in biosample info
    biosample_handle = Entrez.efetch(db="biosample", id=biosample, redmode="xml") 
    biosample_data = biosample_handle.read()
    biosample_handle.close()
    # parse the xml tree, getting taxonomic info
    root = ET.fromstring(biosample_data)
    biosample_species = root.find(".//OrganismName").text
    if record.annotations['organism'] != biosample_species: # if there's a mismatch
        nuc_list.append(record.name)
        nuc_tax_list.append(record.annotations['organism'])
        biosample_list.append(biosample)
        biosample_tax_list.append(biosample_species)
        mismatches += 1
    total += 1

# write to CSV
rows = zip(nuc_list, nuc_tax_list, biosample_list, biosample_tax_list)
with open("taxon_mismatches.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(["GenBank", "GenBank_species", "BioSample", "BioSample_species"])
    for row in rows:
        writer.writerow(row)

print(f"Found {mismatches} mismatches between GenBank and BioSample taxonomy")
print("Conflicting taxonomy written to taxon_mismatches.csv")
print("")
print(f"Found {no_biosample} GenBank records with no associated BioSample")
print("See GenBank_records_no_BioSample.txt for list of accession numbers")

fetch_handle.close()
