#!/usr/bin/env Rscript
# Copyright 2023 Alessandro Gerada alessandro.gerada@liverpool.ac.uk
library(optparse)
library(tidyverse)
library(parallel)
library(glue)
library(Biostrings)
library(Rcpp)

# directories
supported_modality_filters <- c("all", "mic", "disc")
lockBinding("supported_modality_filters", globalenv())

option_list <- list(
  make_option(c("-f", "--filter"),
              default = "MIC",
              help = "Which genomes to download,
              based on available meta-data. One of: all, MIC or disc. 
              default [MIC]"),
  make_option(c("-d", "--database"),
              default = "PATRIC_genomes_AMR.txt",
              help = "Database of genomes downloaded from PATRIC website: 
              https://www.bv-brc.org/docs/quick_references/ftp.html
              default [script_dir/PATRIC_genomes_AMR.txt]"),
  make_option(c("-o", "--output_directory"),
              default = "data/genomes/patric/",
              help = "Folder to save downloaded genomes. 
              default [data/genomes/patric]"),
  make_option(c("-n", "--n_genomes"),
              default = 0,
              help = "Number of genomes to try to download, where: 
              0 = all, 
              -1 = none, 
              default [0]")
)

args <- parse_args(OptionParser(usage = "%script [options] taxonomic_name",
                                option_list = option_list),
                                positional_arguments = 1)

opt <- args$options

# argument clean-up
opt$filter <- tolower(opt$filter)
opt$filter <- ifelse(opt$filter == "disk", "disc", opt$filter)
taxonomic_name <- args$args[[1]]

if (!opt$filter %in% supported_modality_filters) {
  stop(glue("Unable to recognise --filter {opt$filter}, please use one of: 
  {glue_collapse(supported_modality_filters, sep=', ')}"))
}

# make list of unique genome_id for download
patric_amr_list <- read_delim(opt$database, delim = "\t",
                              col_types = cols(.default = "c"))

# select the appropriate samples based on --filter

filtered_data <- patric_amr_list %>%
  filter(str_detect(genome_name, taxonomic_name))

filtered_data <- filtered_data %>% filter(case_when(
  opt$filter == "mic" & measurement_unit == "mg/L" ~ TRUE,
  opt$filter == "disc" & laboratory_typing_method == "Disk diffusion" ~ TRUE,
  opt$filter == "all" ~ TRUE
))

genome_ids <- unique(filtered_data$genome_id)


if (opt$n_genomes < 0) {
  n_downloads <- 0
} else if (opt$n_genomes > 0 & opt$n_genomes < length(genome_ids)) {
  n_downloads <- opt$n_genomes
} else {
  n_downloads <- length(genome_ids)
}

genome_paths <- glue(
  "ftp://ftp.patricbrc.org/genomes/{genome_ids}/{genome_ids}.fna"
  )

if (!dir.exists(opt$output_directory)) dir.create(opt$output_directory)

i <- 1
failures <- 0
while (i <= n_downloads) {
  target_path <- file.path(opt$output_directory,
                          glue("{genome_ids[[i]]}.fna"))
  if (file.exists(target_path)) {
    print(glue("Genome {genome_paths[[i]]} already exists"))
  } else {
    print(glue("Downloading file {i} of {n_downloads}"))
    tryCatch(download.file(genome_paths[[i]],
                  destfile = target_path,
                  mode = "wb"),
             error = function(e) {
              failures <- failures + 1
              print(glue("Unable to download {genome_ids[[i]]}"))
             }
             )
  }
  i <- i + 1
}

print("Downloads complete")
if (failures > 0) {
  print(glue("Failed to download {failures} out of {n_downloads}"))
}
