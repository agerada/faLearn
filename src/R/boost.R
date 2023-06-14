library(tidyverse)
library(AMR)
library(multidplyr)
library(xgboost)
library(optparse)

option_list <- list(
  make_option(c("-v", "--verbose", action = "store_true", default = FALSE,
                help = "Print extra output")),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help =
                "When >1, indicates number of parallel cores to use [default 1]"),
  make_option(c("-r", "--restrict_antibiotics", type = "character", 
                default = NULL, 
                help =
                  "comma-serated list of antibiotics to restrict processing
                of meta-data. This reduces computational time but may miss out 
                on inherent resitant and inferred phenotypes")),
  make_option(c("-t", "--type", type = "character",
                default = "binary",
                help =
                  "type of model to train, only binary supported currently [default binary]")), 
  make_option(c("-i", "--iterations", type = "integer",
                default = 100,
                help =
                  "number of iterations of xgboost model [default 100]")),
  make_option(c("-s", "--save", type = "character", 
                default = FALSE, 
                help =
                  "path to store model file"))
)

args <- parse_args(
  OptionParser(usage = "%script [options] data_file meta_data antibiotic_target", 
               option_list = option_list,
               positional_arguments = 3)
)

opt <- args$options
input_data_path <- args$args[[1]]
meta_data_path <- args$args[[2]]
antibiotic_target <- args$args[[3]]

parallel <- ifelse(opt$cores > 1, TRUE, FALSE)

if (!is.null(opt$restrict_antibiotics)) {
  antibiotic_process_restrict <- unlist(strsplit(opt$restrict_antibiotics))
} else {
  antibiotic_process_restrict <- NULL
}
model_type <- opt$type

clean_raw_mic <- function(mic) {
  # remove leading equals sign
  mic <- stringr::str_remove(mic, pattern = "^=+")
  
  # keep only trimethoprim component
  stringr::str_remove(mic, "/(.)+")
}

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

make_meta_data <- function(patric_meta_data, 
                           kmers_database, 
                           measurement_unit_filter = "mg/L") {
  # meta_data = dataframe of PATRIC meta data
  # kmers = dataframe of kmer reads, needs to include genome_id and 
  # kmer reads in wide format, columns must be named kmer_1 .. kmer_n
  # where n = number of kmers read per organism
  
  # Remove isolates where there are multiple readings for the same organism and
  # antibiotic combo. Unable to determine which to use
  database_filtered <- patric_meta_data %>% 
    group_by(genome_id, antibiotic) %>% 
    filter(n() == 1) 
  
  meta_data <- tibble(genome_id = kmers$genome_id)
  meta_data <- inner_join(database_filtered, meta_data) %>%
    filter(measurement_unit == measurement_unit_filter) 
  
  meta_data <- meta_data %>% 
    pivot_wider(id_cols = genome_id, 
                names_from = antibiotic, 
                values_from = measurement) %>% 
    left_join(distinct(patric_meta_data, genome_id, .keep_all = T))
  
  meta_data <- meta_data %>% 
    mutate(across(any_of(antibiotic_names), clean_raw_mic))
}

meta_data <- make_meta_data(database, kmers)
meta_data <- meta_data %>% select(any_of(c('genome_id', 'genome_name', antibiotic_names)))

# process MICs
clean_up_mics <- function(data, abx_names, parallel = TRUE) {
  if(parallel) {
    cores <- parallel::detectCores()
    cluster <- multidplyr::new_cluster(cores)
    cluster_copy(cluster, c("abx_names", "clean_raw_mic"))
    system.time(
      meta_data <- meta_data %>% 
        partition(cluster) %>% 
        mutate(across(any_of(abx_names), clean_raw_mic)) %>% 
        mutate(across(any_of(abx_names), AMR::as.mic)) %>% 
        collect()
    )
  } else {
    system.time(
      meta_data <- meta_data %>% 
        mutate(across(all_of(abx_names))) %>% 
        mutate(across(all_of(abx_names), AMR::as.mic))
    )
  }
  return(meta_data)
}

convert_mic_to_sir <- function(data, abx_names, parallel = TRUE) {
  if(parallel) {
    cores <- parallel::detectCores()
    cluster <- multidplyr::new_cluster(cores)
    cluster_copy(cluster, c("abx_names"))
    
    system.time(
      meta_data_binary <- data %>% 
        partition(cluster) %>% 
        mutate(organism_name = AMR::as.mo(genome_name)) %>% 
        mutate(across(any_of(abx_names), AMR::as.sir)) %>% 
        collect()
    )
  } else {
    system.time(
      meta_data_binary <- data %>% 
        mutate(organism_name = AMR::as.mo(genome_name)) %>% 
        mutate(across(any_of(abx_names), AMR::as.sir))
    )
  }
  return(meta_data_binary)
}

meta_data_clean_mics <- clean_up_mics(meta_data, antibiotic_names)

if(model_type == "binary") {
  meta_data_sir <- convert_mic_to_sir(meta_data_clean_mics, antibiotic_names)
} else {
  message("Only binary model type is implemented")
  quit("ask")
}

combine_kmer_metadata_and_split <- function(meta_data, 
                                            kmer_data, 
                                            antibiotic,
                                            train_frac = 0.8){
  meta_data_antibiotic <- meta_data %>% 
    left_join(kmer_data, by = "genome_id") %>% 
    filter(!is.na({{antibiotic}}))
  
  train_test_dataset <- meta_data_antibiotic %>% 
    mutate(label = if_else({{antibiotic}} == "R", 1, 0)) 
  
  train_test_dataset <- ungroup(train_test_dataset)
  
  train_ds <- train_test_dataset %>% 
    slice_sample(prop = train_frac, replace = FALSE)
  
  test_ds <- anti_join(train_test_dataset, train_ds)
  
  train_ds <- xgb.DMatrix(data = as.matrix(select(train_ds, starts_with("kmer"))), 
                        label = train_ds$label)
  test_ds <- xgb.DMatrix(data = as.matrix(select(test_ds, starts_with("kmer"))), 
                       label = test_ds$label)
  
  return(list('train' = train_ds, 
              'test' = test_ds))
}

combined_data <- combine_kmer_metadata_and_split(meta_data_sir,
                                                 kmers, 
                                                 antibiotic_target)
training_data <- combined_data[['train']]
testing_data <- combined_data[['test']]

# start with just CIP
make_xgmodel_binary <- function(training_data, iterations = 100, 
                         train_frac = 0.8, cores){ 
  bstDMatrix <- xgboost(data = training_data, max.depth = 2, eta = 1, nthread = cores, 
                        nrounds = iterations, objective = "binary:logistic")
}


xgmodel <- make_xgmodel_binary(training_data,
                               iterations = opt$iterations,
                               cores = opt$cores)

if (opt$save) {
  save(xgmodel, file = opt$save)
}

accuracy_binary <- function(xg_model, test_dataset){
  labels <- getinfo(test_dataset, 'label')
  return(sum(as.numeric(predict(xg_model, test_dataset) > 0.5) == labels) / nrow(test_dataset))
}

accuracy_binary(xgmodel, testing_data)
