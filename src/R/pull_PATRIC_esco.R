library(tidyverse)
library(parallel)
library(glue)
library(Biostrings)
library(Rcpp)
library(bench)

# make list of unique genome_id for download
patric_amr_list <- read_delim('public_data/PATRIC_genomes_AMR.txt', delim='\t')

# select only samples with broth or agar dilution
patric_amr_list_esco_agar_broth <- patric_amr_list %>% 
  filter(str_detect(genome_name, "Escherichia coli")) %>% 
  filter(measurement_unit == "mg/L")

esco_genome_ids <- unique(patric_amr_list_esco_agar_broth$genome_id)

old_wd <- getwd()
setwd("public_data/genomes")

n_downloads <- length(esco_genome_ids)
# check we are not trying to downloading more than available
n_downloads <- ifelse(length(esco_genome_ids) < n_downloads, 
                      length(esco_genome_ids), n_downloads)

esco_genome_paths <- glue("ftp://ftp.patricbrc.org/genomes/{esco_genome_ids}/{esco_genome_ids}.fna")

i <- 1
while(i <= n_downloads){
  print(glue("Downloading file {i} of {n_downloads}"))
  tryCatch(download.file(esco_genome_paths[[i]], 
                destfile = glue("{esco_genome_ids[[i]]}.fna"), 
                mode="wb"), 
           error = function(e) print(glue("Unable to download {esco_genome_ids[[i]]}"))
           )
  i <- i+1
}
setwd(old_wd)



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
