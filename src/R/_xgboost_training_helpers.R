require(tidyverse)
require(AMR)
require(multidplyr)

clean_raw_mic <- function(mic) {
  # remove leading equals sign
  mic <- stringr::str_remove(mic, pattern = "^=+")
  
  # keep only trimethoprim component
  stringr::str_remove(mic, "/(.)+")
}

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

clean_up_mics <- function(data, abx_names, cores = 1) {
  if(cores > 1) {
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

convert_mic_to_sir <- function(data, abx_names, cores = 1) {
  if(cores > 1) {
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

sir_to_binary <- function(data, cols) {
  data %>% mutate(
    across(any_of(cols), function(x) if_else(x == "R", 1, 0)))
}

combine_kmer_meta_data <- function(meta_data,
                                   kmer_data,
                                   antibiotic_filter = NA,
                                   joining_col = "genome_id") {
  meta_data <- meta_data %>%
    left_join(kmer_data, by = joining_col)

  if (antibiotic_filter %in% colnames(meta_data)) {
    antibiotic_filter <- sym(antibiotic_filter)
    meta_data <- meta_data %>% filter(!is.na({{antibiotic_filter}}))
  }
  return(meta_data)
}

split_combined <- function(meta_kmer_data, antibiotic, train_frac = 0.8) {
  train_test_dataset <- ungroup(meta_kmer_data)
  train_ds <- train_test_dataset %>% 
    slice_sample(prop = train_frac, replace = FALSE)
  train_ds_names <- train_ds$genome_id
  train_ds_labels <- train_ds[[antibiotic]]
  

  test_ds <- anti_join(train_test_dataset, train_ds)
  test_ds_names <- test_ds$genome_id
  test_ds_labels <- test_ds[[antibiotic]]
  train_ds <- xgb.DMatrix(data = as.matrix(select(train_ds, starts_with("kmer"))), 
                          label = train_ds_labels)
  test_ds <- xgb.DMatrix(data = as.matrix(select(test_ds, starts_with("kmer"))), 
                         label = test_ds_labels)
  
  return(list('train' = train_ds,
              'train_names' = train_ds_names, 
              'test' = test_ds, 
              'test_names' = test_ds_names))
}

make_xgmodel_binary <- function(training_data, iterations = 100, cores = 1) {
  bstDMatrix <- xgboost(
    data = training_data,
    max.depth = 2,
    eta = 1,
    nthread = cores,
    nrounds = iterations,
    objective = "binary:logistic")
}

make_xgmodel_regression <- function(
  training_data,
  iterations = 100,
  cores = 1
  ) {
    bstDMatrix <- xgboost(
      data = training_data,
      max.depth = 2,
      eta = 1,
      nthread = cores,
      nrounds = iterations,
      objective = "reg:squarederror"
    )
  }

accuracy_binary <- function(xg_model, test_dataset, genome_names = NA){
  labels <- getinfo(test_dataset, 'label')
  summary_table <- tibble(genome = genome_names,
                          label = labels)
  sig_out <- predict(xg_model, test_dataset)
  summary_table$sigmoid_prediction <- sig_out
  predictions <- ifelse(sig_out > 1, 1, 0)
  summary_table$prediction <- predictions
  return(summary_table)
  return(sum(as.numeric(predict(xg_model, test_dataset) > 0.5) == labels) / nrow(test_dataset))
}

mic_double_and_halve <- function(mic) {
  # case_when generates many warnings here because it has to evaluate
  # all RHS conditions for vectorisation
  case_when(
    str_detect(mic, ">") ~ as.mic(mic * 2),
    str_detect(mic, "<") ~ as.mic(mic / 2),
    .default = mic)
}

mic_categorical_to_numeric <- function(data, cols) {
  data <- data %>%
    ungroup %>% 
    mutate(across(any_of(cols), mic_double_and_halve)) %>%
    mutate(across(any_of(cols), as.numeric))
  return(data)
}