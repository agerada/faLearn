strip_filename <- function(paths, extension = "fna") {
  return(
    paths |>
      stringr::str_remove(paste0("\\.", extension, ".*$")) |>
      stringr::str_remove("^.+\\/")
  )
}

batch_filename <- function(dir, batch_num, suffix, ext, batch_lead = 4) {
  # make filename for saving batches
  # batch_num = integer
  # suffix = string descriptor of run, e.g., "k11_meta_data" indicates meta data for
  # 11-kmer run
  # batch_lead = including leading zeros, num of "chars" of batch_num, e.g.,
  # for batch_lead = 4 would result 0001, 0002, etc.
  batch_num_with_lead <- sprintf(paste0("%0", batch_lead, "d"), batch_num)
  path <- file.path(dir, paste0(batch_num_with_lead, suffix, ext))
  return(path)
}

kmers_to_dataframe <- function(kmer_data, k) {
  ids <- names(kmer_data)
  kmer_string <- kmer_data[[1]]$kmer_string
  kmer_data <- data.table::setDT(lapply(kmer_data, \(x) x$kmer_value))
  kmer_data <- data.table::transpose(kmer_data)
  old_names <- names(kmer_data)
  data.table::setnames(kmer_data, old_names, kmer_string)
  kmer_data <- kmer_data[, "genome_id" := ids]
  data.table::setcolorder(kmer_data, "genome_id")
  return(kmer_data)
}

write_sparse_batch <- function(x, file_path, labels = NULL, overwrite = TRUE) {
  if (!is.null(labels) && length(labels) != length(x)) {
    stop("Length of labels and features must be equal")
  }
  if (file.exists(file_path)) {
    if (overwrite) {
      message(paste("Overwriting", file_path))
      file.remove(file_path)
    } else {
      warning(paste(file_path, "already exists. Appending to file."))
    }
  }
  file_conn <- file(file_path, open = "a")
  if (is.null(labels)) {
    labels <- rep(0, length(x))
  }
  for (i in seq_along(x)) {
    features <- paste0(x[[i]]$kmer_index, ":", x[[i]]$kmer_value)
    features <- paste(features, collapse = " ")
    line <- paste(labels[[i]], features)
    writeLines(line, file_conn)
  }
  close(file_conn)
}

write_non_sparse_batch <- function(
    data, target_path, k, overwrite = FALSE,
    file_format = "csv") {

  writing_function <- switch(
    file_format,
    csv = function(x, path) {
      write.csv(x, path, row.names = FALSE)
    },
    parquet = function(x, path) {
      arrow::write_parquet(x, path)
    }
  )
  if (file.exists(target_path) && !overwrite) {
    message("Skipping backup as file is already present.
    Force overwrite if required")
    return()
  }

  if (file.exists(target_path)) {
    message(paste("Overwriting kmer backup to", target_path))
    writing_function(data, target_path)
    return()
  }

  message(paste("Saving kmer backup to", target_path))
  writing_function(data, target_path)
}

genome_paths_from_dir <- function(path,
                                  n_genomes = NULL,
                                  ext = "fna",
                                  random_shuffle = FALSE) {
  confirmed_genomes_paths <- list.files(path,
                                        pattern = paste0("\\.", ext, "$"),
                                        full.names = TRUE)

  if (random_shuffle) {
    message("Shuffling genomes..")
    confirmed_genomes_paths <- sample(confirmed_genomes_paths)
  }

  if (length(confirmed_genomes_paths) < 1) {
    stop("No appropriate files found in input dir")
  }

  n_genomes_to_process <- ifelse(is.null(n_genomes),
                                 length(confirmed_genomes_paths),
                                 n_genomes)

  confirmed_genomes_paths <- confirmed_genomes_paths[1:n_genomes_to_process]

  if (n_genomes_to_process > length(confirmed_genomes_paths)) {
    stop("Number of genomes to process greater than genomes available")
  }

  return(confirmed_genomes_paths)
}

split_paths <- function(paths, split) {
  split_indices <- parallel::splitIndices(
    length(paths), split)

  confirmed_genomes_paths_split <- lapply(split_indices, function(x) {
    sapply(x, function(y) paths[[y]])
  })
  return(confirmed_genomes_paths_split)
}
