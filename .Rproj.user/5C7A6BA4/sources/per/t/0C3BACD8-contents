library(tidyverse)
library(AMR)
library(multidplyr)
library(xgboost)
library(optparse)

source("src/R/_xgboost_training_helpers.R")

option_list <- list(
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = FALSE,
              help = "Print extra output"),
  make_option(c("-c", "--cores"), 
              type = "integer", 
              default = 1, 
              help =
                "When >1, indicates number of parallel cores to use [default 1]"),
  make_option(c("-r", "--restrict_antibiotics"), 
              type = "character",
              default = NULL,
              help =
              "comma-serated list of antibiotics to restrict processing
              of meta-data. This reduces computational time but may miss out 
              on inherent resitant and inferred phenotypes"),
  make_option(c("-t", "--type"), 
              type = "character",
              default = "binary",
              help =
                "type of model to train, only binary supported currently [default binary]"), 
  make_option(c("-i", "--iterations"),
              type = "integer",
              default = 100,
              help =
                "number of iterations of xgboost model [default 100]"),
  make_option(c("-s", "--save"),
              type = "character",
              help =
                "path to store model file")
)

args <- parse_args(
  OptionParser(usage = "%script [options] data_file meta_data antibiotic_target", 
               option_list = option_list),
               positional_arguments = 3
)

opt <- args$options
input_data_path <- args$args[[1]]
meta_data_path <- args$args[[2]]
antibiotic_target <- args$args[[3]]

if (!is.null(opt$restrict_antibiotics)) {
  antibiotic_process_restrict <- unlist(strsplit(opt$restrict_antibiotics,
                                                split = ","))
} else {
  antibiotic_process_restrict <- NULL
}
model_type <- opt$type

database_path <- meta_data_path
kmers_path <- input_data_path
database <- read_delim(database_path, col_types = cols(.default = "c"))

if(!is.null(antibiotic_process_restrict)){
  antibiotic_names <- antibiotic_process_restrict
} else {
  antibiotic_names <- unique(database$antibiotic)
}

kmers <- read_csv(kmers_path,
                  col_types = cols(V1 = col_character())) 

names(kmers) <- c("genome_id", paste0("kmer_", seq(ncol(kmers) - 1)))
kmers <- distinct(kmers, genome_id, .keep_all = TRUE)

meta_data <- make_meta_data(database, kmers)
meta_data <- meta_data %>% select(any_of(c('genome_id', 'genome_name', antibiotic_names)))

# process MICs
meta_data_clean_mics <- clean_up_mics(meta_data, antibiotic_names, 
                                      cores = opt$cores)

if(model_type == "binary") {
  meta_data_sir <- convert_mic_to_sir(meta_data_clean_mics, antibiotic_names,
                                      cores = opt$cores)
} else {
  message("Only binary model type is implemented")
  quit("ask")
}

combined_data <- combine_kmer_metadata_and_split(meta_data_sir,
                                                 kmers, 
                                                 antibiotic_target)
training_data <- combined_data[['train']]
testing_data <- combined_data[['test']]

xgmodel <- make_xgmodel_binary(training_data,
                               iterations = opt$iterations,
                               cores = opt$cores)

if (!is.null(opt$save)) {
  save(xgmodel, file = opt$save)
}

print(accuracy_binary(xgmodel, testing_data))
