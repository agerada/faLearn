library(tidyverse)
library(xgboost)
library(multidplyr)
library(AMR)

source("src/R/_xgboost_training_helpers.R")
source("src/R/_data_converters.R")

database_path <- "/volatile/agerada/molecularMIC/db_backups/esco_all.csv"
#kmers_path <- "/volatile/agerada/molecularMIC/kmers/e_coli_mic/11/rds"
kmers_path <- "/volatile/agerada/molecularMIC/sims/esco_25922/kmers/11"
libsvm_out_path <- "/volatile/agerada/molecularMIC/kmers/e_coli_mic/11/libsvm_gent_mic/"
model_save_path <- "/volatile/agerada/molecularMIC/models/e_coli/11_e_coli_mic_gent_log2.model"
abx <- "gentamicin"
model_type <- "mic"
kmers_format <- "rds"
train_test_split <- 0.8
cores <- 56
train <- FALSE
process_files <- FALSE
patric_db <- load_pre_processed_patric_db(database_path)

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

get_phenotype_value <- switch(model_type, "mic" = "mic")

if (process_files) {
    for (i in seq_along(kmers_filepaths)) {
        batch_name <- strip_filename(kmers_filepaths[[i]], "rds")
        message(paste("Reading ", batch_name))
        x <- readRDS(kmers_filepaths[[i]])
        ids <- names(x)
        phenos <- sapply(
            ids,
            function(x) {
                get_phenotype(
                    x, patric_db, abx, value = get_phenotype_value)
            },
            simplify = FALSE) 
        x <- subset(x, !is.na(phenos))
        phenos <- subset(phenos, !is.na(phenos))

        phenos <- phenos |>
            unlist() |>
            AMR::as.mic() |>
            mic_double_and_halve() |>
            as.numeric()

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

}

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

test_genome_ids <- list_filenames(file.path(libsvm_out_path, "test_meta"), input_format = "csv")
test_genome_ids <- sort(test_genome_ids)
test_genome_ids <- lapply(test_genome_ids, function(x) read_csv(x, col_types = cols(.default = "c")))
test_genome_ids <- unlist(sapply(test_genome_ids, function(x) x$genome_id))



if (train) {
    message("loading training data")
    train_data <- xgb.DMatrix(train_path)
    log_transformed_train_labels <- log2(getinfo(train_data, "label"))
    setinfo(train_data, "label", log_transformed_train_labels)

    message("training model")
    model <- xgboost(
        data = train_data,
        nrounds = 100,
        nthread = cores,
        max.depth = 6,
        eta = 0.3
    )
    xgb.save(model, model_save_path)
} else {
    model <- xgb.load(model_save_path)
}

test_data <- xgb.DMatrix(test_path)
test_labels <- getinfo(test_data, "label")

predictions <- 2 ^ (predict(model, test_data))

results <- tibble(true_data = test_labels, predictions = predictions)
results <- results |>
    mutate(mic_error_ratio = true_data / predictions) |>
    mutate(mic_error_category = case_when(
        mic_error_ratio > 4 ~ "error",
        mic_error_ratio < 0.5 ~ "error",
        TRUE ~ "concordant"
    )) |>
    mutate(true_binary = AMR::as.sir(AMR::as.mic(true_data), mo = as.mo("E. coli"), ab = abx))
