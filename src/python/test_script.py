#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename:      test_script.py
# Author:        Alessandro Gerada
# Date:          28/11/2022
# Copyright:     Alessandro Gerada 2022
# Email:         alessandro.gerada@liverpool.ac.uk

"""
Documentation
"""

import pandas as pd
from src.python.kmer import generate_kmer_perms
esco2018sens = pd.read_csv("../../data/moradigaravand_2018/esco_sens.csv")
esco2018sens = esco2018sens.rename(columns={"Lane.accession": "run_accession"})
esco2018ftp = pd.read_csv("../../data/moradigaravand_2018/esco_ftp.txt", delimiter="\t")

esco2018joined = esco2018sens.merge(esco2018ftp, how="left", on="run_accession")

print(len(generate_kmer_perms(12)))
#[ print(i) for i in (generate_kmer_perms(3))]
"""
esco_seqs = []
with open('data/SAMEA2204418.contigs.fa.gz.filtered.fa') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        esco_seqs.append(record.seq)
print(sum([len(x) for x in esco_seqs]))
esco_kmers = {}
for seqs in esco_seqs:
    esco_kmers = read_kmers(seqs, 10, input_kmer_dict=esco_kmers)
print(esco_kmers)
"""