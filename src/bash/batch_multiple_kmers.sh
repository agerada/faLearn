#!/usr/bin/env bash
# arguments = kmers to batch
for i in $@
do 
    echo "working on $i kmers"
    Rscript ../R/batch_kmer_generator.R -c 4 -k $i -n 10 ../../data/genomes/patric ../../temp
done