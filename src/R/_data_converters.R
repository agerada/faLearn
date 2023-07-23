backup_kmer_data <- function(
  data, target_path, overwrite = FALSE,
  file_format = "csv") {
  # Data must be in the following format:
  # dataframe (or equivalent) with cols
  # genome_id, V1 .. Vn (for n kmers)
  # supports csv or parquet format

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
  for (i in seq_along(x)) {
    features <- paste0(x[[i]]$kmer_index, ":", x[[i]]$kmer_value)
    features <- paste(features, collapse = " ")
    line <- paste(labels[[i]], features)
    writeLines(line, file_conn)
  }
  close(file_conn)
}

load_pre_processed_patric_db <- function(path) {
  return(
    readr::read_csv(path, col_types = cols(.default = "c")) |>
    dplyr::mutate(across(any_of(
      c("resistant_phenotype", "mic_sir", "disk_sir")), AMR::as.sir)) |>
    dplyr::mutate(amr_org = AMR::as.mo(amr_org)) |>
    dplyr::mutate(mic = AMR::as.mic(mic)) |>
    dplyr::mutate(zone = AMR::as.disk(zone)) |>
    # unique genome_id abx combos, prefer MIC
    dplyr::arrange(measurement_unit) |>
    dplyr::group_by(genome_id, antibiotic) |>
    dplyr::slice_head(n = 1)
  )
}

list_filenames <- function(path, input_format, as_kmer_paths = TRUE) {
  # as_kmer_paths returns kmer_paths object

  files_list <- list.files(
    path,
    pattern = paste0("\\.", input_format, "$"),
    full.names = TRUE,
    ignore.case = TRUE)
  files_list <- normalizePath(files_list) # remove double slashes
  if (as_kmer_paths) {
    class(files_list) <- "kmer_paths"
  }
  return(files_list)
}

strip_filename <- function(paths, extension) {
  return(
    paths |>
      stringr::str_remove(paste0("\\.", extension, ".*$")) |>
      stringr::str_remove("^.+\\/")
  )
}

sort.kmer_paths <- function(x, decreasing = FALSE) {
  extension <- sub(".*\\.(.*?)$", "\\1", x)
  filename_prefix <- as.integer(
    stringr::str_extract(strip_filename(x, extension), "\\d+$"))
  sorting_indices <- sort(
    filename_prefix, decreasing = decreasing, index.return = TRUE)$ix

  x <- x[sorting_indices]
  class(x) <- "kmer_paths"
  return(x)
}

train_test_filesystem <- function(
  path_to_files,
  file_ext, split = 0.8, train_folder = "train", test_folder = "test") {
  libsvm_filepaths <- list_filenames(path_to_files, file_ext)

  if (
    length(libsvm_filepaths) == 0 &&
    dir.exists(file.path(path_to_files, train_folder)) &&
    dir.exists(file.path(path_to_files, test_folder))) {
      message("files already appear to be in train test subdirectories")
      out_paths <- normalizePath(file.path(path_to_files, c(train_folder, test_folder)))
      names(out_paths) <- c(train_folder, test_folder)
      return(out_paths)
    }

  libsvm_filepaths <- sort(libsvm_filepaths)

  dir.create(file.path(path_to_files, train_folder))
  dir.create(file.path(path_to_files, test_folder))

  splitting_index <- split * length(libsvm_filepaths)
  train_libsvm_paths <- head(libsvm_filepaths, splitting_index)
  test_libsvm_paths <- tail(libsvm_filepaths, length(libsvm_filepaths) - splitting_index)

  target_ext <- paste0(".", file_ext)

  file.rename(
      from = train_libsvm_paths,
      to = file.path(
          path_to_files,
          train_folder,
          paste0(strip_filename(train_libsvm_paths, file_ext), target_ext)))

  file.rename(
      from = test_libsvm_paths,
      to = file.path(
          path_to_files,
          test_folder,
          paste0(strip_filename(test_libsvm_paths, file_ext), target_ext)))

  out_paths <- normalizePath(file.path(path_to_files, c(train_folder, test_folder)))
  names(out_paths) <- c(train_folder, test_folder)
  return(out_paths)
}