#!/usr/bin/env bash
IN_PATH="/volatile/agerada/molecularMIC/sims/esco_25922/mutations"
OUT_PATH="/volatile/agerada/molecularMIC/sims/esco_25922/mutations/annotations"

for i in "$IN_PATH"/*.fna
do 
    file_name=$(basename -s .fna ${i})
    amrfinder -n $i -O Escherichia -o ${IN_PATH}${file_name}.tsv
done