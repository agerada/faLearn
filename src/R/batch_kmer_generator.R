#!/usr/bin/env Rscript
# Copyright 2022 Alessandro Gerada alessandro.gerada@liverpool.ac.uk
library(optparse)

option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output"), 
  make_option(c("-c", "--cores"), type="integer", default=1, 
              help="When >1, indicates number of parallel cores to use [default 1]"), 
  make_option(c("-k", "--kmers"), type="integer", default=3, 
              help="Kmer length [default 3]"), 
  make_option(c("-n","--n_genomes"), type="integer", 
              help="Max number of genomes to process [default all]"), 
  make_option(c("-p", "--split"), type="integer", default=1, 
              help="Split data into s chunks for parallel processing [default 1]"), 
  make_option(c("-a", "--anchor"), action="store_true", default=FALSE, 
              help="Include unobserved permutations"), 
  make_option(c("-s", "--simplify"), action="store_true", default=FALSE,
              help="Store only the kmer counts (key: value -> value)")
)
print(getwd())
args <- parse_args(OptionParser(usage = "%script [options] input_dir output_dir", 
                               option_list=option_list), 
                  positional_arguments=2)

opt <- args$options
dirs <- args$args
input_dir <- dirs[[1]]
output_dir <- dirs[[2]]

confirmed_genomes_paths <- list.files(input_dir, pattern="\\.fna$", full.names = TRUE)
confirmed_genomes_ids <- confirmed_genomes_paths |> stringr::str_remove(".fna.*$")

if (length(confirmed_genomes_paths) < 1) { stop("No .fna files found in input dir")}

n_genomes_to_process <- ifelse(is.null(opt$n_genomes), 
                               length(confirmed_genomes_paths), 
                               opt$n_genomes)

confirmed_genomes_paths <- confirmed_genomes_paths[1:n_genomes_to_process]
confirmed_genomes_ids <- confirmed_genomes_ids[1:n_genomes_to_process]

if (n_genomes_to_process > length(confirmed_genomes_paths)) stop("Number of genomes to process 
                                                          greater than genomes available")

path_to_rcpp_script <- file.path( dirname(getwd()), "Rcpp/kmer.cpp")
Rcpp::sourceCpp(path_to_rcpp_script)

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

convert_to_kmers <- function(x) {
  if(opt$verbose) print(paste("Working on",x))
  x <- as.character(unlist(Biostrings::readDNAStringSet(x)))
  x <- kmers_pointed(x, kmer=opt$kmers, anchor=opt$anchor, simplify=opt$simplify)
  x
}

split_indices <- parallel::splitIndices(length(confirmed_genomes_paths), opt$split)
confirmed_genomes_paths_split <- lapply(split_indices, function(x){
  sapply(x, function(y) confirmed_genomes_paths[[y]])
})
confirmed_genomes_ids_split <- lapply(split_indices, function(x){
  sapply(x, function(y) confirmed_genomes_ids[[y]])
})

if (opt$cores > 1) {
  if (Sys.info()['sysname'] != 'Windows') {
    for (i in seq_along(confirmed_genomes_paths_split))
    {
      if (opt$verbose) message(paste("Working on chunk", i, "of", opt$split))
      output <- parallel::mclapply(confirmed_genomes_paths_split[[i]], 
                         convert_to_kmers, mc.cores = opt$cores)
      names(output) <- confirmed_genomes_ids_split[[i]]
      save(output, file=file.path(output_dir, 
                                  paste0(opt$kmers, "_kmer_data", i, ".RData")))
      remove(output)
    }
  }
  else{
    # if running on Windows, we require a few extra steps for parallel to work
    # use snow package instead
    cl <- snow::makeCluster(opt$cores, type='SOCK')
    # next expose required variables and Rcpp scripts. Note that any package 
    # function calls should be done without namespace (therefore e.g. Rcpp::)
    # otherwise run library(package) on all cores
    snow::clusterExport(cl, 'opt$kmers')
    snow::clusterEvalQ(cl, Rcpp::sourceCpp(path_to_rcpp_script))
    for (i in seq_along(confirmed_genomes_paths_split))
    {
      if (opt$verbose) message(paste("Working on chunk", i, "of", opt$split))
      output <- pbapply::pblapply(confirmed_genomes_paths_split[[i]], 
                         convert_to_kmers, cl=cl)
      names(output) <- confirmed_genomes_ids_split[[i]]
      save(output, file=file.path(output_dir, 
                                  paste0(opt$kmers, "_kmer_data", i, ".RData")))
      remove(output)
    }
    snow::stopCluster(cl)
  }
} else {
  for (i in seq_along(confirmed_genomes_paths_split))
  {
    if (opt$verbose) message(paste("Working on chunk", i, "of", opt$split))
    output <- lapply(confirmed_genomes_paths_split[[i]], convert_to_kmers)
    names(output) <- confirmed_genomes_ids_split[[i]]
    save(output, file=file.path(output_dir, 
                                paste0(opt$kmers, "_kmer_data", i, ".RData")))
    remove(output)
  }
}

