library(tidyverse)
library(xgboost)
library(multidplyr)
library(AMR)

source("src/R/_xgboost_training_helpers.R")

# options
opt <- list()
opt$restrict_antibiotics <- "ciprofloxacin,gentamicin,meropenem"
opt$type <- "regression"
opt$cores <- 56
opt$iterations <- 100

model_type <- opt$type
if (!is.null(opt$restrict_antibiotics)) {
  antibiotic_process_restrict <- unlist(strsplit(opt$restrict_antibiotics,
                                                split = ","))
} else {
  antibiotic_process_restrict <- NULL
}

# paths
input_data_path <- "/volatile/agerada/molecularMIC/kmers/acinetobacter/6/6_kmer_backup.csv"
meta_data_path <- "/home/agerada/molecularMIC/data/databases/patric/PATRIC_genomes_AMR.txt"
antibiotic_target <- "ciprofloxacin"


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

if (model_type == "binary") {
  meta_data_sir <- convert_mic_to_sir(meta_data_clean_mics, antibiotic_names,
                                      cores = opt$cores)
  meta_data_binary <- sir_to_binary(meta_data_sir, antibiotic_names)
  meta_data <- meta_data_binary
} else if (model_type == "regression") {
  meta_data_mic_numeric <- mic_categorical_to_numeric(meta_data_clean_mics,
                                                      antibiotic_names)
  meta_data <- meta_data_mic_numeric
}

combined_data <- combine_kmer_meta_data(meta_data, kmers, antibiotic_target)
split_data <- split_combined(combined_data, "ciprofloxacin")

training_data <- split_data[['train']]
testing_data <- split_data[['test']]

if (model_type == "binary") {
  xgmodel <- make_xgmodel_binary(training_data,
                                iterations = opt$iterations,
                                cores = opt$cores)
} else if (model_type == "regression") {
  xgmodel <- make_xgmodel_regression(
    training_data,
    iterations = 1000,
    cores = opt$cores
  )
}

if (!is.null(opt$save)) {
  save(xgmodel, file = opt$save)
}
