strip_filename <- function(paths) {
  tools::file_path_sans_ext(basename(paths))
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
#'
#' @examples
#' set.seed(123)
#' # create 10 random DNA files
#' tmp_dir <- tempdir()
#' # remove any existing .fna files
#' file.remove(
#'   list.files(tmp_dir, pattern = "*.fna", full.names = TRUE)
#' )
#'
#' for (i in 1:10) {
#'  writeLines(paste0(">", i, "\n", paste0(sample(c("A", "T", "C", "G"),
#'  100, replace = TRUE), collapse = "")), file.path(tmp_dir, paste0(i, ".fna")))
#' }
#'
#' # split files into train and test directories
#' paths <- train_test_filesystem(tmp_dir,
#'                                file_ext = "fna",
#'                                split = 0.8,
#'                                shuffle = TRUE,
#'                                overwrite = TRUE)
#'
#' list.files(paths[["train"]])
#' list.files(paths[["test"]])
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

    if (!overwrite) {
      out_paths <- normalizePath(file.path(path_to_files, c(train_folder, test_folder)))
      names(out_paths) <- c(train_folder, test_folder)
      return(out_paths)
    }

  }

  if (!isTRUE(shuffle)) {
    libsvm_filepaths <- sort(libsvm_filepaths)
  }

  if (dir.exists(file.path(path_to_files, train_folder)) | dir.exists(file.path(path_to_files, test_folder))) {
    message("train or test folders already exist")
    if (!overwrite) {
      stop("Aborting, force using overwrite = TRUE")
    }
  }

  unlink(file.path(path_to_files, train_folder), recursive = TRUE)
  unlink(file.path(path_to_files, test_folder), recursive = TRUE)

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
#' @description
#' This function reorganises files that have been split into train and test
#' directories using train_test_filesystem() back into a single directory.
#' This is a convenience function to reverse the effects of
#' train_test_filesystem().
#'
#' @export
#'
#' @return Logical vector, indicated success or failure for each file
#'
#' @examples
#' set.seed(123)
#' # create 10 random DNA files
#' tmp_dir <- tempdir()
#' # remove any existing .fna files
#' file.remove(
#'  list.files(tmp_dir, pattern = "*.fna", full.names = TRUE)
#' )
#' for (i in 1:10) {
#' writeLines(paste0(">", i, "\n", paste0(sample(c("A", "T", "C", "G"),
#'   100, replace = TRUE), collapse = "")), file.path(tmp_dir, paste0(i, ".fna")))
#' }
#'
#' # split files into train and test directories
#' paths <- train_test_filesystem(tmp_dir,
#'                                file_ext = "fna",
#'                                split = 0.8,
#'                                shuffle = TRUE,
#'                                overwrite = TRUE)
#' # combine files back into a single directory
#' combined_file_system(tmp_dir, "fna")
#' list.files(tmp_dir)
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

  train_outcomes <- file.rename(
    from = train_files,
    to = target_train_paths)

  test_outcomes <- file.rename(
    from = test_files,
    to = target_test_paths)

  return(all(c(train_outcomes, test_outcomes)))
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
#'
#' @return Logical vector, indicating success or failure for each file
#'
#' @description
#' This is simply a wrapper around file.copy/file.rename that allows for
#' filtering by a logical vector (move_which). This can replicate the behaviour
#' of a predicate function (see example), and may be easier to read.
#'
#' @examples
#' set.seed(123)
#' # create 10 random DNA files
#' tmp_dir <- tempdir()
#' # remove any existing .fna files
#' file.remove(
#'  list.files(tmp_dir, pattern = "*.fna", full.names = TRUE)
#' )
#' for (i in 1:10) {
#' writeLines(paste0(">", i, "\n", paste0(sample(c("A", "T", "C", "G"),
#'  100, replace = TRUE), collapse = "")), file.path(tmp_dir, paste0(i, ".fna")))
#' }
#'
#' # move files with even numbers to a new directory
#' new_dir <- file.path(tempdir(), "even_files")
#' unlink(new_dir, recursive = TRUE)
#' move_files(tmp_dir,
#'            new_dir,
#'            move_which = as.integer(
#'               tools::file_path_sans_ext(
#'                   list.files(tmp_dir, pattern = "*.fna"))) %% 2 == 0,
#'            ext = "fna")
#' list.files(new_dir)
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
    return(
      sapply(filtered_paths, \(x) {
        outcome <- file.copy(from = x,
                  to = file.path(target_dir, basename(x)))
        p()
        outcome
      })
    )
  }
  return(
    sapply(filtered_paths, \(x) {
      outcome <- file.rename(from = x,
                  to = file.path(target_dir, basename(x)))
      p()
      outcome
    })
  )
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
#' @param path_to_files path containing files or vector of filepaths
#' @param file_ext file extension to filter
#' @param split train-test split
#' @param train_target_path name of train file to save as (by default, will be
#' train.txt in the path_to_files directory)
#' @param test_target_path name of test file to save as (by default, will be
#' test.txt in the path_to_files directory)
#' @param names_backup name of file to save backup of filename metadata (by
#' default, will be names.csv in the path_to_files directory)
#' @param shuffle randomise prior to splitting
#' @param overwrite overwrite target files
#'
#' @return named list of paths to created train/test files, original filenames
#'
#' @description
#' This function combines files into a train and test set, stored on disk.
#' It can be used in combination with genomes_to_kmer_libsvm() to create a
#' dataset that can be loaded into XGBoost (either by first creating an
#' xgboost::DMatrix, or by using the data argument in xgboost::xgb.train()
#' or xgboost::xgb.cv()). The following three files will be created:
#'
#' 1) train.txt - the training data
#' 2) test.txt - the testing data (if split < 1)
#' 3) names.csv - a csv file containing the original filenames and their
#'   corresponding type (train or test)
#'
#' The function will check if the data is already in the appropriate format
#' and will not overwrite unless forced using the overwrite argument.
#'
#' By providing 1.0 to the split argument, the function can be used to combine
#' files without a train-test split. In this case, all the files will be
#' classed as 'train', and there will be no 'test' data. This is useful if
#' one wants to perform cross-validation using xgboost::xgb.cv() or
#' MIC::xgb.cv.lowmem(). It is also possible to combine all data into train
#' and then perform splitting after loading into an xgboost::DMatrix, using
#' xgboost::slice().
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' # create 10 random libsvm files
#' tmp_dir <- tempdir()
#' # remove any existing .txt files
#' file.remove(
#' list.files(tmp_dir, pattern = "*.txt", full.names = TRUE)
#' )
#' for (i in 1:10) {
#'  # each line is K: V
#'  writeLines(paste0(i, ": ", paste0(sample(1:100, 10, replace = TRUE),
#'  collapse = " ")), file.path(tmp_dir, paste0(i, ".txt")))
#'  }
#'
#'  # split files into train and test directories
#'  paths <- split_and_combine_files(
#'   tmp_dir,
#'   file_ext = "txt",
#'   split = 0.8,
#'   train_target_path = file.path(tmp_dir, "train.txt"),
#'   test_target_path = file.path(tmp_dir, "test.txt"),
#'   names_backup = file.path(tmp_dir, "names.csv"),
#'   overwrite = TRUE)
#'
#'  readLines(paths[["train"]])
split_and_combine_files <- function(path_to_files,
                                    file_ext = ".txt",
                                    split = 0.8,
                                    train_target_path = NULL,
                                    test_target_path = NULL,
                                    names_backup = NULL,
                                    shuffle = TRUE,
                                    overwrite = FALSE) {
  file_ext <- gsub("^\\.", "", file_ext)

  train_target_path <- ifelse(is.null(train_target_path),
                              file.path(path_to_files, "train.txt"),
                              train_target_path)
  test_target_path <- ifelse(is.null(test_target_path),
                             file.path(path_to_files, "test.txt"),
                             test_target_path)
  names_backup <- ifelse(is.null(names_backup),
                         file.path(path_to_files, "names.csv"),
                         names_backup)

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

  p <- progressr::progressor(along = c(train_libsvm_paths, test_libsvm_paths))
  for (file in train_libsvm_paths) {
    content <- readLines(file, warn = FALSE)
    write(content, train_target_path, append = TRUE, sep = "/n")
    p(glue::glue("Processing {file}"))
  }

  for (file in test_libsvm_paths) {
    content <- readLines(file, warn = FALSE)
    write(content, test_target_path, append = TRUE, sep = "/n")
    p(glue::glue("Processing {file}"))
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
