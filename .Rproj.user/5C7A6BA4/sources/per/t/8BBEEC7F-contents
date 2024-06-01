library(tidyverse)
library(xgboost)
library(multidplyr)
library(AMR)

source("src/R/_xgboost_training_helpers.R")
source("src/R/_data_converters.R")

kmers_path <- "/volatile/agerada/molecularMIC/sims/esco_25922/kmers/11/libsvm"
meta_data_path <- "/volatile/agerada/molecularMIC/sims/esco_25922/kmers/11/meta_data"
model_load_path <- "/volatile/agerada/molecularMIC/models/e_coli/11_e_coli_mic_gent_log2.model"
kmers_format <- "txt"
meta_data_format <- "csv"

kmers_filepaths <- list_filenames(kmers_path, kmers_format, as_kmer_paths = F)
kmers_filepaths <- sort(kmers_filepaths)

kmers <- xgb.DMatrix(kmers_path)

meta_data_filepaths <- list_filenames(meta_data_path, meta_data_format, as_kmer_paths = F)

model <- xgb.load(model_load_path)

predictions <- 2 ^ (predict(model, kmers))
genome_ids <- read_csv(meta_data_filepaths)$genome_id
names(predictions) <- genome_ids
