library(tidyverse)
library(parallel)
library(glue)
library(Biostrings)
library(Rcpp)
library(bench)
sourceCpp("src/Rcpp/kmer.cpp", rebuild=T)

# options
parallel <- FALSE
kmer_var <- 10
n_genomes_to_process <- 4

# directories
database_txt <- "data/databases/patric/PATRIC_genomes_AMR.txt"
genomes_dir <- "data/genomes/patric"
kmer_output_dir <- "data/output/kmers"

# get metadata for genomes
patric_amr_metadata <- read_delim(database_txt, delim='\t', 
                                  col_types = cols(.default="c"))
available_genomes_ids <- list.files(genomes_dir) %>% str_remove(".fna.*$")
patric_amr_metadata <- patric_amr_metadata %>% 
  filter(genome_id %in% available_genomes_ids)

# check that we do not have some discrepancy
if(length(base::unique(patric_amr_metadata$genome_id)) != 
   length(available_genomes_ids) ) warning("One or more genomes do not have 
                                           matching metadata..")

# generate list of paths to available genomes with matching metadata
confirmed_genomes_paths <- paste0(patric_amr_metadata$genome_id,".fna")
confirmed_genomes_paths <- confirmed_genomes_paths[confirmed_genomes_paths %in% list.files(genomes_dir)]
confirmed_genomes_paths <- base::unique(confirmed_genomes_paths)
confirmed_genomes_ids <- confirmed_genomes_paths %>% str_remove(".fna.*$")
confirmed_genomes_paths <- file.path(genomes_dir, confirmed_genomes_paths)

if(!dir.exists(kmer_output_dir)) dir.create(kmer_output_dir, recursive = TRUE)

if (parallel) {
  output <- mclapply(confirmed_genomes_paths[1:n_genomes_to_process], 
           function(x) x %>% readDNAStringSet() %>% unlist() %>% 
             as.character() %>% kmers_pointed(., kmer=kmer_var, anchor=T, simplify=T) %>% 
             unlist) 
} else {
  output <- lapply(confirmed_genomes_paths[1:n_genomes_to_process], 
                   function(x) x %>% readDNAStringSet() %>% unlist() %>% 
                     as.character() %>% kmers_pointed(., kmer=kmer_var, anchor=T, simplify=T) %>% 
                     unlist) 
}
names(output) <- confirmed_genomes_ids[1:n_genomes_to_process]
save(output, file=file.path(kmer_output_dir, paste0(kmer_var, "_kmer_data.RData")))
