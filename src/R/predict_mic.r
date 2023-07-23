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

kmers_filepaths <- list_filenames(kmers_path, kmers_format)
kmers_filepaths <- sort(kmers_filepaths)

kmers <- lapply(kmers_filepaths, xgb.DMatrix)

meta_data_filepaths <- list_filenames(meta_data_path, meta_data_format)

model <- xgb.load(model_load_path)

predictions <- lapply(kmers, function(x) predict(model, x))

predictions <- lapply(predictions, function(x) 2 ^ x)
predictions <- unlist(predictions)
genome_ids <- unlist(sapply(meta_data_filepaths, function(x) {
    csv <- read.csv(x)
    csv$genome_id
}))
names(predictions) <- genome_ids
