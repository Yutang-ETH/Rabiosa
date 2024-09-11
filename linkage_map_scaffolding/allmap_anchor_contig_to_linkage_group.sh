#!/bin/bash

# Prepare ALLMAPS input
python -m jcvi.assembly.allmaps merge RSM.csv RSFM.csv -o JM-2.bed -w weights.txt

# Run ALLMAPS
python -m jcvi.assembly.allmaps path JM-2.bed scaffolds_FINAL.fasta -w weights.txt
