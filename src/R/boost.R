library(tidyverse)
library(AMR)
library(multidplyr)
library(xgboost)


parallel <- TRUE

clean_raw_mic <- function(mic) {
  # remove leading equals sign
  mic <- stringr::str_remove(mic, pattern = "^=+")
  
  # keep only trimethoprim component
  stringr::str_remove(mic, "/(.)+")
}

database_path <- "data/databases/patric/PATRIC_genomes_AMR.txt"
kmers_path <- "data/output/projections/6_kmer_backup.csv"

database <- read_delim(database_path, col_types = cols(.default = "c"))
kmers <- read_csv(kmers_path,
                  col_types = cols(V1 = col_character())) 

names(kmers) <- c("genome_id", paste0("kmer_", seq(ncol(kmers) - 1)))
kmers <- distinct(kmers, genome_id, .keep_all = TRUE)

# Remove isolates where there are multiple readings for the same organism and
# antibiotic combo. Unable to determine which to use
database_filtered <- database %>% 
  group_by(genome_id, antibiotic) %>% 
  filter(n() == 1) 

meta_data <- tibble(genome_id = kmers$genome_id)
meta_data <- inner_join(database_filtered, meta_data) %>%
  filter(measurement_unit == "mg/L") 

antibiotic_names <- unique(meta_data$antibiotic)
meta_data <- meta_data %>% 
  pivot_wider(id_cols = genome_id, 
              names_from = antibiotic, 
              values_from = measurement) %>% 
  left_join(distinct(database, genome_id, .keep_all = T))

meta_data <- meta_data %>% 
  mutate(across(all_of(antibiotic_names), clean_raw_mic))

if(parallel) {
  cores <- parallel::detectCores()
  cluster <- multidplyr::new_cluster(cores)
  cluster_copy(cluster, c("antibiotic_names", "clean_raw_mic"))
  system.time(
    meta_data <- meta_data %>% 
      partition(cluster) %>% 
      mutate(across(all_of(antibiotic_names), clean_raw_mic)) %>% 
      mutate(across(all_of(antibiotic_names), AMR::as.mic)) %>% 
      collect()
  )
} else {
  system.time(
    meta_data <- meta_data %>% 
      mutate(across(all_of(antibiotic_names), AMR::as.mic))
  )
}

antibiotic_subset <- c("ciprofloxacin", "gentamicin", "ampicillin")

if(parallel) {
  cores <- parallel::detectCores()
  cluster <- multidplyr::new_cluster(cores)
  cluster_copy(cluster, c("antibiotic_subset"))
  
  system.time(
    meta_data_binary <- meta_data %>% 
      partition(cluster) %>% 
      mutate(organism_name = AMR::as.mo(genome_name)) %>% 
      mutate(across(all_of(antibiotic_subset), AMR::as.sir)) %>% 
      collect()
  )
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

combined_data <- combine_kmer_metadata_and_split(meta_data_binary,
                                                 kmers, 
                                                 ciprofloxacin)
training_data <- combined_data[['train']]
testing_data <- combined_data[['test']]

# start with just CIP
make_xgmodel_binary <- function(training_data, iterations = 100, 
                         train_frac = 0.8){ 
  bstDMatrix <- xgboost(data = training_data, max.depth = 2, eta = 1, nthread = cores, 
                        nrounds = iterations, objective = "binary:logistic")
}
  
xgmodel <- make_xgmodel_binary(training_data)

accuracy_binary <- function(xg_model, test_dataset){
  labels <- getinfo(test_dataset, 'label')
  return(sum(as.numeric(predict(xg_model, test_dataset) > 0.5) == labels) / nrow(test_ds))
}

accuracy_binary(xgmodel, testing_data)
