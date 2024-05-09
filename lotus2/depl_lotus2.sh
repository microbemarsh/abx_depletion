#!/bin/bash

# Load conda env
source /Users/austinmarshall/mambaforge/etc/profile.d/conda.sh
conda activate lotus2

#Trying with conda

lotus2 -amplicon_type SSU -refDB SLV -taxAligner 0 -CL vsearch -t 12 -i /path/to/sequences/ -m /path/to/DEPL_lotus2_map_v2.txt -o depl_lotus2
