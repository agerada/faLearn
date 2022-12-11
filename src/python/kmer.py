#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Filename:      kmer.py
# Author:        Alessandro Gerada
# Date:          28/11/2022
# Copyright:     Alessandro Gerada 2022
# Email:         alessandro.gerada@liverpool.ac.uk

"""
Documentation
"""

def read_kmers(seq, k, input_kmer_dict):
    kmer_dict = input_kmer_dict if input_kmer_dict else {}
    n = len(seq) + 1 - k
    for i in range(n):
        kmer_i = seq[i:i+k]
        if kmer_i in kmer_dict:
            kmer_dict[kmer_i] += 1
        else:
            kmer_dict[kmer_i] = 1
    return kmer_dict
from itertools import product
def generate_kmer_perms(k, bases = ["A", "T", "C", "G"]):
    """
    adapted from https://stackoverflow.com/questions/30340361/permutations-with-repetition-using-recursion-javascript
    There are n**k permutations possible, where n = bases
    """
    output = []
    def recursive_inner(permutation):
        if len(permutation) > k - 1:
            output.append(permutation)
            return
        for i in bases:
            recursive_inner(permutation + i)
    recursive_inner("")
    return output


    """
    return product(bases, repeat=k)
    result = []
    for base in bases:
        if k == 1:
            result.append(base)
        else:
            perms = generate_kmer_perms(k - 1, bases)
            for j in perms:
                result.append(perms)
    return result
    """
