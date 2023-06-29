library(tidyverse)
library(xgboost)
library(multidplyr)
library(AMR)
library(arrow)

source("src/R/_xgboost_training_helpers.R")

# options
opt <- list()
opt$restrict_antibiotics <- "ciprofloxacin,gentamicin,meropenem"
opt$type <- "parquet"
opt$cores <- 56
opt$iterations <- 100
opt$verbose <- TRUE

model_type <- opt$type

if (!is.null(opt$restrict_antibiotics)) {
  antibiotic_process_restrict <- unlist(strsplit(opt$restrict_antibiotics,
                                                split = ","))
} else {
  antibiotic_process_restrict <- NULL
}

# paths
#input_data_path <- "/volatile/agerada/molecularMIC/kmers/acinetobacter/10/split/10_kmer_data1.csv"
input_data_path <- "/volatile/agerada/molecularMIC/kmers/e_coli_all/3/"
meta_data_path <- "/home/agerada/molecularMIC/data/databases/patric/PATRIC_genomes_AMR.txt"
antibiotic_target <- "ciprofloxacin"
input_format <- "csv"

database_path <- meta_data_path
kmers_path <- input_data_path

if (opt$verbose) message("Reading meta database")
database_raw <- read_delim(database_path, col_types = cols(.default = "c"))

if (!is.null(antibiotic_process_restrict)) {
  if (opt$verbose) message(paste("Restricting further antibiotic processing to",antibiotic_process_restrict))
  antibiotic_names <- antibiotic_process_restrict
} else {
  antibiotic_names <- unique(database$antibiotic)
}

if (opt$verbose) message("Reading kmers csv")
kmer_files <- list.files(
  kmers_path,
  pattern = paste0("*.", input_format),
  full.names = TRUE)
input_schema <- schema(genome_id = string())
kmers <- arrow::open_dataset(
  kmer_files,
  format = input_format)
kmers <- data.frame(kmers)
names(kmers) <- c("genome_id", paste0("feature_", seq(ncol(kmers) - 1)))

if (opt$verbose) message("Keeping only unique genome_ids")
kmers <- distinct(kmers, genome_id, .keep_all = TRUE)

if (opt$verbose) message("Creating meta data for kmers")
database <- filter(database_raw, genome_id %in% kmers$genome_id)
if (!is.null(opt$restrict_antibiotics)) {
  database <- filter(database, antibiotic %in% antibiotic_process_restrict)
}
database <- mutate(database, amr_org = AMR::as.mo(genome_name))
database <- mutate(
  database, 
  measurement = if_else(measurement_unit == "mg/L", clean_raw_mic(measurement), measurement))
database <- mutate(database, mic = if_else(measurement_unit == "mg/L", AMR::as.mic(measurement), NA))
database <- mutate(database, zone = if_else(measurement_unit == "mm", AMR::as.disk(measurement), NA))

if (model_type == "binary") {
  database <- populate_sir(database, cores = opt$cores)
  database <- database %>% mutate(resistant_phenotype =
    case_when(
      !is.na(mic_sir) ~ mic_sir,
      !is.na(disk_sir) ~ disk_sir,
      .default = resistant_phenotype
    )
  )
}

#write_csv(database, "/volatile/agerada/molecularMIC/binary_database_bak.csv")

meta_data <- make_meta_data(
  database, kmers, binary = if_else(
    model_type == "binary", TRUE, FALSE))
meta_data <- meta_data %>% select(any_of(c('genome_id', 'genome_name', antibiotic_names)))

if (model_type == "binary") {
  if (opt$verbose) message("Converting SIR to binary")
  meta_data <- sir_to_binary(meta_data, antibiotic_names)
} else if (model_type == "regression") {
  if (opt$verbose) message("Converting MIC to numeric (doubling if >=, halving if <=)")
  meta_data <- mic_categorical_to_numeric(meta_data,
                                                      antibiotic_names)
}

if (opt$verbose) message("Combining meta data and kmers")

system.time(
  combined_data <- combine_kmer_meta_data(meta_data, kmers, antibiotic_target)
)

system.time(
  split_data <- split_combined(
    combined_data, "ciprofloxacin")
)

training_data <- split_data[['train_data']]
testing_data <- split_data[['test_data']]

if (model_type == "binary") {
  xgmodel <- make_xgmodel_binary(training_data,
                                iterations = opt$iterations,
                                cores = opt$cores)
} else if (model_type == "regression") {
  xgmodel <- make_xgmodel_regression(
    training_data,
    iterations = opt$iterations,
    cores = opt$cores)
  }

if (!is.null(opt$save)) {
  save(xgmodel, file = opt$save)
}

predictions <- predict(xgmodel, testing_data)
if (model_type == "binary") {
  predictions <- if_else(predictions < 0.5, 0, 1)
} else if (model_type == "regression") {
  predictions <- if_else(predictions < 0, 0, predictions)
}
plot(predictions, split_data$test_labels)
table(data.frame(predict = predictions, labels = split_data$test_labels))
