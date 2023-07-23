library(tidyverse)
library(xgboost)
library(multidplyr)
library(AMR)

source("src/R/_xgboost_training_helpers.R")
source("src/R/_data_converters.R")

database_path <- "/volatile/agerada/molecularMIC/db_backups/esco_all.csv"
kmers_path <- "/volatile/agerada/molecularMIC/kmers/e_coli_mic/11/rds"
libsvm_out_path <- "/volatile/agerada/molecularMIC/kmers/e_coli_mic/11/test/"
model_save_path <- "/volatile/agerada/molecularMIC/models/e_coli/11_e_coli_mic_cip.model"
abx <- "ciprofloxacin"
kmers_format <- "rds"
train_test_split <- 0.8
cores <- 56
train <- TRUE

patric_db <- load_pre_processed_patric_db(database_path)
patric_db <- sir_to_binary(patric_db, "resistant_phenotype")

kmers_filepaths <- list_filenames(kmers_path, kmers_format)
kmers_filepaths <- sort(kmers_filepaths)

get_phenotype <- function(id, database, ab, value = "resistant_phenotype") {
    pheno_index <- database$genome_id == id & database$antibiotic == ab
    if (!any(pheno_index)) return(NA)
    if (sum(pheno_index) > 1) {
        stop("Multiple phenotypes for genome_id, don't know which to choose.")
    }
    return(database[[value]][pheno_index])
}

if (!dir.exists(libsvm_out_path)) dir.create(libsvm_out_path, recursive = TRUE)

for (i in seq_along(kmers_filepaths)) {
    batch_name <- strip_filename(kmers_filepaths[[i]], "rds")
    message(paste("Reading ", batch_name))
    x <- readRDS(kmers_filepaths[[i]])
    ids <- names(x)
    phenos <- sapply(
        ids,
        function(x) get_phenotype(x, patric_db, abx), simplify = FALSE)
    x <- subset(x, !is.na(phenos))
    phenos <- subset(phenos, !is.na(phenos))

    meta_data <- data.frame(genome_id = names(x))

    meta_data_path <- file.path(
        libsvm_out_path,
        paste0(stringr::str_replace(
            batch_name, "data", "meta_data"), ".csv"))

    message(paste("Writing meta data to", meta_data_path))
    write.csv(meta_data, meta_data_path, row.names = FALSE)

    message(paste("Writing", paste0(batch_name, ".libsvm")))
    write_sparse_batch(x, file.path(libsvm_out_path, paste0(batch_name, ".libsvm")), phenos)
}


if (!train) q()

paths <- train_test_filesystem(libsvm_out_path, "libsvm", train_test_split)

train_path <- paths[["train"]]
test_path <- paths[["test"]]

genome_meta_data <- list_filenames(libsvm_out_path, "csv")
genome_meta_data <- sort(genome_meta_data)

paths <- train_test_filesystem(
    libsvm_out_path, "csv", train_test_split,
    train_folder = "train_meta", test_folder = "test_meta")

train_genome_ids <- list_filenames(file.path(libsvm_out_path, "train_meta"), input_format = "csv")
train_genome_ids <- sort(train_genome_ids)
train_genome_ids <- lapply(train_genome_ids, function(x) read_csv(x, col_types = cols(.default = "c")))
train_genome_ids <- unlist(sapply(train_genome_ids, function(x) x$genome_id))

train_data <- xgb.DMatrix(train_path)

model <- xgboost(
    data = train_data,
    objective = "binary:logistic",
    nrounds = 100,
    nthread = cores,
    max.depth = 6,
    eta = 0.3
)

test_DMatrix <- xgb.DMatrix(test_path)
test_labels <- getinfo(test_DMatrix, "label")

predictions <- predict(model, test_DMatrix)
predictions_binary <- ifelse(predictions > 0.5, 1, 0)

error_rate <- (1 / length(test_labels)) * (sum(test_labels != predictions_binary))

xgb.save(model, model_save_path)