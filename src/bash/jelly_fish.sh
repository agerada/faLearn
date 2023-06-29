#!/usr/bin/env bash
IN_PATH="/volatile/agerada/molecularMIC/genomes/e_coli_all/"
OUT_PATH="/volatile/agerada/molecularMIC/kmers/e_coli_all/3/"

for i in "$IN_PATH"/*.fna
do 
    echo $i
    jellyfish 
done