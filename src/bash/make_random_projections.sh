#!/usr/bin/env bash
# arguments = kmers to batch
for i in $@
do 
    echo "working on $i kmers"
    Rscript ../R/random_projections.R -m -b $i 10 ../../data/output/kmers/${i}_kmer_data1.RData ../../temp
done
