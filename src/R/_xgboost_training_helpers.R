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

combine_kmer_metadata_and_split <- function(meta_data, 
                                            kmer_data, 
                                            antibiotic,
                                            train_frac = 0.8){
  antibiotic <- sym(antibiotic)
  meta_data_antibiotic <- meta_data %>% 
    left_join(kmer_data, by = "genome_id") %>% 
    filter(!is.na({{antibiotic}}))
  
  train_test_dataset <- meta_data_antibiotic %>% 
    mutate(label = if_else({{antibiotic}} == "R", 1, 0)) 
  
  train_test_dataset <- ungroup(train_test_dataset)
  
  train_ds <- train_test_dataset %>% 
    slice_sample(prop = train_frac, replace = FALSE)
  train_ds_names <- train_ds$genome_id
  train_ds <- xgb.DMatrix(data = as.matrix(select(train_ds, starts_with("kmer"))), 
                          label = train_ds$label)
  
  test_ds <- anti_join(train_test_dataset, train_ds)
  test_ds_names <- test_ds$genome_id
  test_ds <- xgb.DMatrix(data = as.matrix(select(test_ds, starts_with("kmer"))), 
                         label = test_ds$label)
  
  return(list('train' = train_ds,
              'train_names' = train_ds_names, 
              'test' = test_ds, 
              'test_names' = test_ds_names))
}

make_xgmodel_binary <- function(training_data, iterations = 100, 
                                train_frac = 0.8, cores){ 
  bstDMatrix <- xgboost(data = training_data, max.depth = 2, eta = 1, nthread = cores, 
                        nrounds = iterations, objective = "binary:logistic")
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
