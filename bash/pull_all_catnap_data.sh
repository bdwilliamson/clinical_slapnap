#!/bin/bash

# pull all CATNAP data
mkdir -p dat/catnap

# pull CATNAP data from LANL
wget -O dat/catnap/assay.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/assay.txt"
wget -O dat/catnap/viruses.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/viruses.txt"
wget -O dat/catnap/virseqs_aa.fasta "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/virseqs_aa.fasta"
wget -O dat/catnap/abs.txt "https://www.hiv.lanl.gov/cgi-bin/common_code/download.cgi?/scratch/NEUTRALIZATION/abs.txt"
