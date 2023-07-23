library(tidyverse)
library(xgboost)
library(multidplyr)
library(AMR)

source("src/R/_xgboost_training_helpers.R")
source("src/R/_data_converters.R")

kmers_path <- "/volatile/agerada/molecularMIC/sims/esco_25922/kmers/11"
model_load_path <- "/volatile/agerada/molecularMIC/models/e_coli/11_e_coli_mic_gent_log2.model"
kmers_format <- "txt"

kmers_filepaths <- list_filenames(kmers_path, kmers_format)
kmers_filepaths <- sort(kmers_filepaths)

kmers <- xgb.DMatrix(kmers_filepaths)

model <- xgb.load(model_load_path)

predictions <- predict(model, kmers)
predictions <- 2 ^ predictions
