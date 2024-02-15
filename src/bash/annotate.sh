#!/usr/bin/env bash
IN_PATH="/volatile/agerada/molecularMIC/genomes/cgr/"
OUT_PATH="/volatile/agerada/molecularMIC/annotations/cgr/"

for i in "$IN_PATH"/*.fasta
do 
    file_name=$(basename -s .fasta ${i})
    amrfinder -n $i -O Escherichia -o ${OUT_PATH}${file_name}.tsv
done