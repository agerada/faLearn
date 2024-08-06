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
      utils::write.csv(x, path, row.names = FALSE)
    },
    parquet = function(x, path) {
      rlang::check_installed("arrow", "Writing parquet format needs arrow library")
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

#' Organise files into a train-test filesystem
#'
#' @param path_to_files directory containing files
#' @param file_ext file extension to filter
#' @param split training data split
#' @param train_folder name of training folder (subdirectory), will be created
#' if does not exist
#' @param test_folder name of testing folder (subdirectory), will be created
#' if does not exist
#' @param shuffle randomise files when splitting (if FALSE, files will be
#' sorted by filename prior to splitting)
#' @param overwrite force overwrite of files that already exist
#'
#' @return named vector of train and test directories
#' @export
train_test_filesystem <- function(path_to_files,
                                  file_ext,
                                  split = 0.8,
                                  train_folder = "train",
                                  test_folder = "test",
                                  shuffle = TRUE,
                                  overwrite = FALSE) {
  file_ext <- gsub("^\\.", "", file_ext)
  libsvm_filepaths <- list.files(path_to_files,
                                 pattern = paste0("*.", file_ext),
                                 full.names = TRUE,
                                 ignore.case = TRUE)
  if (isTRUE(shuffle)){
    libsvm_filepaths <- sample(libsvm_filepaths)
  }

  if (
    length(libsvm_filepaths) == 0 &&
    dir.exists(file.path(path_to_files, train_folder)) &&
    dir.exists(file.path(path_to_files, test_folder))) {
    message("files already appear to be in train test subdirectories")
    out_paths <- normalizePath(file.path(path_to_files, c(train_folder, test_folder)))
    names(out_paths) <- c(train_folder, test_folder)
    return(out_paths)
  }

  if (!isTRUE(shuffle)) {
    libsvm_filepaths <- sort(libsvm_filepaths)
  }

  dir.create(file.path(path_to_files, train_folder))
  dir.create(file.path(path_to_files, test_folder))

  splitting_index <- split * length(libsvm_filepaths)
  train_libsvm_paths <- utils::head(libsvm_filepaths, splitting_index)
  test_libsvm_paths <- utils::tail(libsvm_filepaths, length(libsvm_filepaths) - splitting_index)

  target_ext <- paste0(".", file_ext)

  target_train_paths <- file.path(
    path_to_files,
    train_folder,
    basename(train_libsvm_paths))

  target_test_paths <- file.path(
    path_to_files,
    test_folder,
    basename(test_libsvm_paths))

  if (any(file.exists(target_train_paths)) | any(file.exists(target_test_paths))) {
    message("The following target paths already exist:")
    sapply(c(target_train_paths[file.exists(target_train_paths)],
             target_test_paths[file.exists(target_test_paths)]),
           message)
    if (!overwrite) {
      stop("Aborting, force using overwrite = TRUE")
    }
  }

  file.rename(
    from = train_libsvm_paths,
    to = target_train_paths)

  file.rename(
    from = test_libsvm_paths,
    to = target_test_paths)

  out_paths <- normalizePath(file.path(path_to_files, c(train_folder, test_folder)))
  names(out_paths) <- c(train_folder, test_folder)
  return(out_paths)
}

#' Combine train and test filesystem into single folder
#'
#' @param path_to_folders path containing test and train folders; files will be
#' moved here
#' @param file_ext file extension to filter
#' @param train_folder train folder subdirectory name
#' @param test_folder test folder subdirectory name
#' @param overwrite force overwrite of files that already exist
#'
#' @export
combined_file_system <- function(path_to_folders,
                                 file_ext,
                                 train_folder = "train",
                                 test_folder = "test",
                                 overwrite = FALSE) {
  file_ext <- gsub("^\\.", "", file_ext)
  train_files <- list.files(file.path(path_to_folders, train_folder),
                            pattern = paste0("*.", file_ext),
                            full.names = TRUE,
                            ignore.case = TRUE)
  test_files <- list.files(file.path(path_to_folders, test_folder),
                           pattern = paste0("*.", file_ext),
                           full.names = TRUE,
                           ignore.case = TRUE)
  if (length(train_folder) == 0) {
    warning("No files found in train folder.")
  }
  if (length(test_folder) == 0) {
    warning("No files found in test folder.")
  }

  target_train_paths <- file.path(
    path_to_folders,
    basename(train_files))

  target_test_paths <- file.path(
    path_to_folders,
    basename(test_files))

  if (any(file.exists(target_train_paths)) | any(file.exists(target_test_paths))) {
    message("The following target paths already exist:")
    sapply(c(target_train_paths[file.exists(target_train_paths)],
             target_test_paths[file.exists(target_test_paths)]),
           message)
    if (!overwrite) {
      stop("Aborting, force using overwrite = TRUE")
    }
  }

  file.rename(
    from = train_files,
    to = target_train_paths)

  file.rename(
    from = test_files,
    to = target_test_paths)
  return(NULL)
}

#' Move or copy files using logical vector
#'
#' @param source_dir move from directory
#' @param target_dir move to directory
#' @param move_which logical vector to filter (or use TRUE to move all)
#' @param ext file extension to filter
#' @param copy copy files (rather than move)
#'
#' @export
move_files <- function(source_dir,
                       target_dir,
                       move_which,
                       ext = ".txt",
                       copy = FALSE) {
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  ext <- gsub("^\\.", "", ext)
  file_paths <- list.files(source_dir,
                           pattern = paste0("*.", ext),
                           full.names = TRUE,
                           ignore.case = TRUE)
  filtered_paths <- subset(file_paths, move_which)
  p <- progressr::progressor(along = filtered_paths)
  if (isTRUE(copy)) {
    sapply(filtered_paths, \(x) {
      file.copy(from = x,
                to = file.path(target_dir, basename(x)))
      p()
    })
    return()
  }
  sapply(filtered_paths, \(x) {
    file.rename(from = x,
                to = file.path(target_dir, basename(x)))
    p()
  })
  return()
}

#' Removes multiple slashes in a path or url
#'
#' @param path character vector
#'
#' @return character vector of paths without duplicate slashes
replace_multiple_slashes <- function(path) {
  url_pattern <- "^(https?://)"

  cleaned_text <- gsub("([^:/])(/{2,})", "\\1/", path)

  is_url <- grepl(url_pattern, path)
  cleaned_text[is_url] <- sub(url_pattern, "\\1", cleaned_text[is_url])

  return(cleaned_text)
}

#' Create test train files from a number of files
#'
#' @param path_to_files path containing files
#' @param file_ext file extension to filter
#' @param split train-test split
#' @param train_target_path name of train file to save as
#' @param test_target_path name of test file to save as
#' @param names_backup name of file to save backup of filename metadata
#' @param shuffle randomise prior to splitting
#' @param overwrite overwrite target files
#'
#' @return named list of paths to created train/test files, original filenames
#' @export
split_and_combine_files <- function(path_to_files,
                                    file_ext = ".txt",
                                    split = 0.8,
                                    train_target_path = file.path(getwd(),
                                                                  "train.txt"),
                                    test_target_path = file.path(getwd(),
                                                                 "test.txt"),
                                    names_backup = file.path(getwd(),
                                                             "names.csv"),
                                    shuffle = TRUE,
                                    overwrite = FALSE) {
  file_ext <- gsub("^\\.", "", file_ext)

  attempted_load <- is_test_train_combined(train_path = train_target_path,
                                           test_path = test_target_path,
                                           names_path = names_backup)
  if (!is.null(attempted_load) & !overwrite) {
    message(
      "Data already seems to be in the appropriate format.
No changes made. Use overwrite to force changes.")
    return(attempted_load)
  }

  all_target_files <- c(train_target_path,
                        test_target_path,
                        names_backup)

  if (any(
    file.exists(all_target_files))
    & !overwrite
  ) stop("Target files already exist, use overwrite to force.")

  suppressWarnings(file.remove(all_target_files))

  if (length(path_to_files) == 1 & all(dir.exists(path_to_files))) {
    libsvm_filepaths <- list.files(path_to_files,
                                   pattern = paste0("*.", file_ext),
                                   full.names = TRUE,
                                   ignore.case = TRUE)
  } else if (is.character(path_to_files)) {
    libsvm_filepaths <- path_to_files
    if (any(!endsWith(libsvm_filepaths, file_ext))) {
      warning(paste("path_to_files contains files that do not have ext:", file_ext))
    }
  } else {
    stop("path_to_files must be directory or character vector of filepaths")
  }

  if (isTRUE(shuffle)){
    libsvm_filepaths <- sample(libsvm_filepaths)
  }

  libsvm_filepaths <- subset(libsvm_filepaths,
                             !(libsvm_filepaths %in% c(train_target_path,
                                                       test_target_path)))

  splitting_index <- split * length(libsvm_filepaths)
  train_libsvm_paths <- utils::head(libsvm_filepaths, splitting_index)
  test_libsvm_paths <- utils::tail(libsvm_filepaths, length(libsvm_filepaths) - splitting_index)

  train_filenames <- basename(train_libsvm_paths)
  test_filenames <- basename(test_libsvm_paths)

  sapply(dirname(all_target_files), function(x) {
    if (!dir.exists(x)) dir.create(x)
  })
  file.create(all_target_files)

  p <- progressr::progressor(along = train_libsvm_paths)
  for (file in train_libsvm_paths) {
    content <- readLines(file, warn = FALSE)
    write(content, train_target_path, append = TRUE, sep = "/n")
    p()
  }

  p <- progressr::progressor(along = test_libsvm_paths)
  for (file in test_libsvm_paths) {
    content <- readLines(file, warn = FALSE)
    write(content, test_target_path, append = TRUE, sep = "/n")
    p()
  }

  readr::write_csv(data.frame(type = c(rep("train", length(train_filenames)),
                                       rep("test", length(test_filenames))),
                              name = c(train_filenames,test_filenames)),
                   names_backup)

  return(list("train" = replace_multiple_slashes(train_target_path),
              "test" = replace_multiple_slashes(test_target_path),
              "train_names" = replace_multiple_slashes(train_filenames),
              "test_names" = replace_multiple_slashes(test_filenames)))
}

is_test_train_combined <- function(train_path, test_path, names_path) {
  if (all(file.exists(c(train_path, test_path, names_path)))) {
    tryCatch({
      names_data <- readr::read_csv(names_path,
                                    col_types = readr::cols(.default = "c"))
      train_files <- names_data$name[names_data$type == "train"]
      test_files <- names_data$name[names_data$type == "test"]
    }, error = function(e) {
      message("Data appears to already be in a combined test-train split, but
              error on loading meta-data:")
      message(conditionMessage(e))
    })
    return(list("train" = replace_multiple_slashes(train_path),
                "test" = replace_multiple_slashes(test_path),
                "train_names" = replace_multiple_slashes(train_files),
                "test_names" = replace_multiple_slashes(test_files)))
  }
  return(NULL)
}
