library(tidyverse)
library(parallel)
library(glue)
library(Biostrings)
library(Rcpp)
library(bench)

seq1 <- readDNAStringSet(filepath = "public_data/genomes/562.101663.fna")
seq2 <- readDNAStringSet(filepath = "public_data/SAMEA2204418.contigs.fa.gz.filtered.fa")

sourceCpp('kmer.cpp')
kmers(s)
s <- as.character(seq1[[1]])
bench::mark(
  kmers(as.character(unlist(seq1)), kmer=10), iterations = 1
)
read_kmer_r(s, kmer=10)

kmers(as.character(unlist(seq1)), kmer=10)
bench::mark(
  kmers_pointed(s, 10, anchor = T), 
  kmers_pointed(s, 10, anchor = F), 
  kmers(s, 10), 
  iterations=10, 
  check=F
)


mclapply(list(as.character(unlist(seq1)),
              as.character(unlist(seq2))), function(x) kmers(x, kmer=10), 
         mc.cores = 4)

# make kmer data
sequences_filenames <- list.files(path="public_data/genomes", pattern = "*.fna", 
                                  full.names = TRUE)
# try 50 sequences first
test_kmers <- mclapply(sequences_filenames[1:4], 
                       function(x) {
                         seq <- Biostrings::readDNAStringSet(x)
                         seq <- as.character(unlist(seq))
                         kmers(seq, kmer=10)
                       })

# try new function
bench::mark(
  mclapply(sequences_filenames[1:4], 
           function(x) {
             seq <- Biostrings::readDNAStringSet(x)
             seq <- as.character(unlist(seq))
             kmers_pointed(seq, kmer=10, anchor=T)[[2]]
           }), 
  memory=FALSE, 
  iterations = 1
)
