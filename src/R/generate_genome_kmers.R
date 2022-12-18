library(tidyverse)
library(parallel)
library(glue)
library(Biostrings)
library(Rcpp)
library(bench)
library(snow)

sourceCpp("src/Rcpp/kmer.cpp", rebuild=T)

# options
parallel <- T
kmer_var <- 3
n_genomes_to_process <- 4
clusters_var <- 8

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

convert_to_kmers <- function(x) {
  print(paste("Working on",x))
  x <- as.character(unlist(Biostrings::readDNAStringSet(x)))
  x <- kmers_pointed(x, kmer=kmer_var, anchor=T, simplify=F)
  x
}

if (parallel) {
  if (Sys.info()['sysname'] != 'Windows') {
  output <- mclapply(confirmed_genomes_paths[1:n_genomes_to_process], 
           convert_to_kmers)
  }
  else{
    # if running on Windows, we require a few extra steps for parallel to work
    # use snow package instead
    cl <- snow::makeCluster(clusters_var, type='SOCK')
    # next expose required variables and Rcpp scripts. Note that any package 
    # function calls should be done without namespace (therefore e.g. Rcpp::)
    # otherwise run library(package) on all cores
    clusterExport(cl, 'kmer_var')
    clusterEvalQ(cl, Rcpp::sourceCpp("src/Rcpp/kmer.cpp"))
    output <- clusterApply(cl, confirmed_genomes_paths[1:n_genomes_to_process], 
                 convert_to_kmers)
    stopCluster(cl)
  }
} else {
  output <- lapply(confirmed_genomes_paths[1:n_genomes_to_process], 
                   function(x) x %>% readDNAStringSet() %>% unlist() %>% 
                     as.character() %>% kmers_pointed(., kmer=kmer_var, anchor=T, simplify=T) %>% 
                     unlist) 
}

names(output) <- confirmed_genomes_ids[1:n_genomes_to_process]
save(output, file=file.path(kmer_output_dir, paste0(kmer_var, "_kmer_data.RData")))
