#!/usr/bin/env Rscript
# Copyright 2022 Alessandro Gerada alessandro.gerada@liverpool.ac.uk
library(optparse)
source(here::here("src/R", "_data_converters.R"))

option_list <- list( 
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print extra output"),
  make_option(c("-c", "--cores"), type = "integer", default = 1,
              help =
              "When >1, indicates number of parallel cores to use [default 1]"),
  make_option(c("-n", "--n_genomes"), type = "integer",
              help = "Max number of genomes to process [default all]"),
  make_option(c("-o", "--overwrite"), action = "store_true", default = FALSE,
              help =
              "Overwrite output file if already present [default false]"),
  make_option(c("-m", "--matrix"), action = "store_true", default = FALSE,
              help =
              "Save matrix (in output_dir) used for 
              random projections [default false]"),
  make_option(c("-b", "--backup"), action = "store_true", default = FALSE,
              help =
              "Make backup (in output_dir) of kmers 
              in csv format [default false]")
  )

args <- parse_args(OptionParser(usage =
  "%script [options] kmers_in dimensions_out input_Rdata_file output_dir",
                                option_list = option_list),
                   positional_arguments = 4)

opt <- args$options
kmers <- as.numeric(args$args[[1]])
dimensions_out <- as.numeric(args$args[[2]])
input_file <- args$args[[3]]
output_dir <- args$args[[4]]

if (!file.exists(input_file)) stop("Input .RData file not found")
print(paste("Input file:", input_file))

if (!grepl(".RData$", input_file)) stop("Input file must be .RData extension")
if (!dir.exists(output_dir)) stop("Output directory does not exist")
data_env <- new.env()

# Load data
object_name <- tryCatch(
  invisible(load(input_file, envir = data_env)),
  error = function(e) {
    message("Error loading .RData input file: ")
    message(e)
    return(NA)
  },
  warning = function(w) {
    message("Attempt to load .RData input file generated warning: ")
    message(w)
  },
  finally = {
    message("Loaded .RData input file")
  }
)
if (length(object_name) > 1) {
  stop("Only one object per .RData file is supported")
}
data_list <- as.list(data_env)
data <- data_list[[object_name]]

data <- lapply(data, function(x) {
  x <- unlist(x)
  names(x) <- NULL
  x})

kmer_perms <- 4 ^ kmers
# remove genomes that include N
data <- Filter(\(x) length(x) == kmer_perms, data)

# set n of genomes to process
n_to_process <- ifelse(is.null(opt$n_genomes), length(data), opt$n_genomes)

message("Generating random projections")
set.seed(25)
kmer_data <- matrix(0, nrow = n_to_process, ncol = kmer_perms)

for (i in seq_along(head(data, n_to_process))) {
  for (j in 1:kmer_perms){
    kmer_data[i, j] <- data[[i]][j]
    }
}

if (opt$backup) {
  backup_file <- file.path(output_dir, paste0(kmers, "_kmer_backup.csv"))
  kmer_Rdata_to_csv(cbind(names(data), kmer_data), backup_file, 
                    overwrite = opt$overwrite)
}

random_matrix <- matrix(
  rnorm(kmer_perms * dimensions_out), nrow = dimensions_out, ncol = kmer_perms)

X <- kmer_data
R <- random_matrix

if (opt$matrix) {
  matrix_backup <- file.path(
    output_dir, paste0(kmers, "_kmer_random_matrix.csv"))
  if (file.exists(matrix_backup) && !opt$overwrite) {
    message("Random matrix backup - file already exists, skipping. 
    Use --overwrite to force overwrite")
  } else if (file.exists(matrix_backup)) {
    message(paste("Overwriting random matrix backup to", matrix_backup))
    write.csv(random_matrix, matrix_backup, row.names = FALSE)
  } else {
    message(paste("Saving random matrix backup to", matrix_backup))
    write.csv(random_matrix, matrix_backup, row.names = FALSE)
  }
}

normalise <- function(x) x / sqrt(sum(x ^ 2))
random_proj <- R %*% t(X)
random_proj <- (1 / sqrt(length(head(data, n_to_process)))) * random_proj

random_proj <- t(random_proj)
random_proj <- cbind(names(data), random_proj)

message("Saving output..")
output_file <- file.path(
  output_dir, paste0(kmers, "_kmer_random_projections.csv"))
if (file.exists(output_file) && !opt$overwrite) {
  stop(paste("Output file already exists:", output_file))
} else if (file.exists(output_file)) {
  message(paste("Overwriting output file to:", output_file))
  write.csv(random_proj, output_file, row.names = FALSE)
  message(paste("Output saved to", output_file))
} else {
  message(paste("Saving output file to:", output_file))
  write.csv(random_proj, output_file, row.names = FALSE)
  message(paste("Output saved to", output_file))
  }
