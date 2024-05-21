max_batch_split <- 9999 # filenames break if greater


verbose_message <- function(s, verbosity = FALSE){
  if (verbosity) message(s)
}

flat_stringset <- function(x) {
  as.character(unlist(Biostrings::readDNAStringSet(x)))
}

genome_to_kmers <- function(x,
                             kmers =3,
                             simplify =T,
                             anchor= T,
                             drop_n=T) {
  message(paste("Working on", x))
  x <- flat_stringset(x)
  kmers(x,
        k = kmers,
        anchor = anchor,
        simplify = simplify,
        clean_up = drop_n,
        key_as_int = FALSE)
}

convert_to_int_indexed_kmers <- function(x,
                                         kmers =3,
                                         simplify =T,
                                         anchor= T,
                                         drop_n=T) {
  message(paste("Working on", x))
  x <- flat_stringset(x)
  kmers(x,
        k = kmers,
        anchor = anchor,
        simplify = simplify,
        clean_up = drop_n,
        key_as_int = TRUE)
}



genomes_to_kmers_libsvm <- function(input_dir,
                                    output_dir,
                                    cores = 1,
                                    kmers = 3,
                                    n_genomes = NULL,
                                    split = 1,
                                    anchor = FALSE,
                                    simplify = FALSE,
                                    drop_n = FALSE,
                                    integer_index = FALSE,
                                    random_shuffle = FALSE) {
  if (split > max_batch_split) {
    stop(paste("Max number of batches =", max_batch_split))
  }

  confirmed_genome_paths <- genome_paths_from_dir(input_dir,
                                                  n_genomes = n_genomes,
                                                  ext = "fna",
                                                  random_shuffle = random_shuffle)
  confirmed_genome_ids <- strip_filename(confirmed_genome_paths)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  split_genome_paths <- split_paths(confirmed_genome_paths, split)
  split_genome_ids <- split_paths(confirmed_genome_ids, split)

  if (cores == 1) {
    for (i in seq_along(split_genome_paths)) {
      message(paste("Working on chunk", i, "of", split))
      output <- lapply(split_genome_paths[[i]], convert_to_int_indexed_kmers,
                       kmers = kmers,
                       simplify = simplify,
                       anchor= anchor,
                       drop_n=drop_n)
      names(output) <- split_genome_ids[[i]]

      save(output, file = batch_filename(
        output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
      saveRDS(output, file = batch_filename(
        output_dir, i, paste0("_k", kmers, "_kmer_data"), ".rds"))

      write_sparse_batch(
        output,
        batch_filename(
          output_dir, i, paste0("_k", kmers, "_libsvm"), ".txt"))

      meta_data <- data.frame(genome_id = names(output))
      write.csv(
        meta_data,
        batch_filename(
          output_dir, i, paste0("_k", kmers, "_meta_data"), ".csv"))
      remove(output)
    }
  }

}
#
# kmer_batch <- function(input_dir,
#                        output_dir,
#                        verbose = FALSE,
#                        cores = 1,
#                        kmers = 3,
#                        n_genomes = NULL,
#                        split = 1,
#                        anchor = FALSE,
#                        simplify = FALSE,
#                        drop_n = FALSE,
#                        formats = "libsvm",
#                        integer_index = FALSE,
#                        random_shuffle = FALSE) {
#   # max_batch_split <- 9999 # filenames break if greater
#   save_formats <- tolower(unlist(strsplit(formats, split = ",")))
# #
# #   if (split > max_batch_split) {
# #     stop(paste("Max number of batches =", max_batch_split))
# #   }
# #
# #   confirmed_genomes_paths <- list.files(input_dir,
# #                                         pattern = "\\.fna$", full.names = TRUE)
# #
# #   if (random_shuffle) {
# #     verbose_message("Shuffling genomes..", verbose)
# #     confirmed_genomes_paths <- sample(confirmed_genomes_paths)
# #   }
# #
# #   confirmed_genomes_ids <- strip_filename(confirmed_genomes_paths, "fna")
# #
# #   if (length(confirmed_genomes_paths) < 1) {
# #     stop("No .fna files found in input dir")
# #   }
#
#   # n_genomes_to_process <- ifelse(is.null(n_genomes),
#   #                                length(confirmed_genomes_paths),
#   #                                n_genomes)
#   #
#   # confirmed_genomes_paths <- confirmed_genomes_paths[1:n_genomes_to_process]
#   # confirmed_genomes_ids <- confirmed_genomes_ids[1:n_genomes_to_process]
# #
# #   if (n_genomes_to_process > length(confirmed_genomes_paths)) {
# #     stop("Number of genomes to process greater than genomes available")
# #   }
#
#   # if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#   #
#   # split_indices <- parallel::splitIndices(
#   #   length(confirmed_genomes_paths), split)
#   #
#   # confirmed_genomes_paths_split <- lapply(split_indices, function(x) {
#   #   sapply(x, function(y) confirmed_genomes_paths[[y]])
#   # })
#   # confirmed_genomes_ids_split <- lapply(split_indices, function(x) {
#   #   sapply(x, function(y) confirmed_genomes_ids[[y]])
#   # })
#
#   if (cores > 1) {
#     if (Sys.info()["sysname"] != "Windows") {
#       for (i in seq_along(confirmed_genomes_paths_split))
#       {
#
#         if ("libsvm" %in% save_formats || integer_index) {
#           verbose_message(
#             paste("Working on int indexed chunk", i, "of", split),
#             verbose)
#           # libsvm format needs integer index insted of kmer string, so different
#           # call to Rcpp library
#           # also needs integer index if -i flagged
#           output <- parallel::mclapply(confirmed_genomes_paths_split[[i]],
#                                        convert_to_int_indexed_kmers,
#                                        mc.cores = cores,
#                                        simplify = simplify)
#           names(output) <- confirmed_genomes_ids_split[[i]]
#           if ("libsvm" %in% save_formats) {
#             # write_sparse_batch(
#             #   output,
#             #   batch_filename(
#             #     output_dir, i, paste0("_k", kmers, "_libsvm"), ".txt"))
#
#             # write meta data separately for libsvm
#           #   meta_data <- data.frame(genome_id = names(output))
#           #   write.csv(
#           #     meta_data,
#           #     batch_filename(
#           #       output_dir, i, paste0("_k", kmers, "_meta_data"), ".csv"))
#           # }
#
#           # if ("rdata" %in% save_formats) {
#           #   save(output, file = batch_filename(
#           #     output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#           # }
#           #
#           # if ("rds" %in% save_formats) {
#           #   rds_filepath <- batch_filename(
#           #     output_dir, i, paste0("_k", kmers, "_kmer_data"), ".rds")
#           #   verbose_message(paste("Writing rds file to", rds_filepath), verbose)
#           #   saveRDS(output, file = rds_filepath)
#           # }
#
#           remove(output)
#
#           if (integer_index) {
#             # formats below are not integer indexed
#             next()
#           }
#         }
#
#         non_integer_indexed_formats <- c("rdata", "rds", "csv", "parquet")
#         if (any(save_formats %in% non_integer_indexed_formats)) {
#           verbose_message(
#             paste("Working on chunk", i, "of", split),
#             verbose)
#
#           output <- parallel::mclapply(confirmed_genomes_paths_split[[i]],
#                                        convert_to_kmers, mc.cores = cores)
#           names(output) <- confirmed_genomes_ids_split[[i]]
#         }
#
#         if ("rdata" %in% save_formats) {
#           save(output, file = batch_filename(
#             output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#         }
#
#         if ("rds" %in% save_formats) {
#           saveRDS(output, file = batch_filename(
#             output_dir, i, paste0("_k", kmers, "_kmer_data"), ".rds"))
#         }
#
#         if ("csv" %in% save_formats || "parquet" %in% save_formats) {
#           # further processing required for these formats
#           if (!simplify) {
#             message("CSV or parquet save without --simplify is not implemented")
#           } else {
#             output <- lapply(output, function(x) {
#               x <- unlist(x)
#               x})
#
#             # keep only items with length == 4 ^ kmers
#             target_length <- 4 ^ kmers
#             output <- subset(
#               output,
#               sapply(output, function(x) {
#                 ifelse(length(x) == target_length, TRUE, FALSE)
#               }))
#
#             genome_ids <- as.character(names(output))
#             data.table::setDT(output)
#             output <- data.table::transpose(output)
#             output <- data.table(genome_id = genome_ids, output)
#             colnames(output) <- c("genome_id", paste0(
#               "kmer_", seq(from = 1, to = ncol(output) - 1)))
#             output[, genome_id := genome_ids]
#
#             if ("parquet" %in% save_formats) {
#               backup_kmer_data(
#                 output,
#                 batch_filename(
#                   output_dir,
#                   i,
#                   paste0("_k", kmers, "_kmer_data"),
#                   ".parquet"),
#                 file_format = "parquet")
#             }
#
#             if ("csv" %in% save_formats) {
#               backup_kmer_data(
#                 output,
#                 batch_filename(
#                   output_dir,
#                   i,
#                   paste0("_k", kmers, "_kmer_data"),
#                   ".csv"),
#                 overwrite = TRUE)
#             }
#
#           }
#           if (exists("output")) remove(output)
#         }
#       }
#     } else {
#       # if running on Windows, we require a few extra steps for parallel to work
#       # use snow package instead
#       cl <- snow::makeCluster(cores, type = 'SOCK')
#       # next expose required variables and Rcpp scripts. Note that any package
#       # function calls should be done without namespace (therefore e.g. Rcpp::)
#       # otherwise run library(package) on all cores
#       snow::clusterExport(cl, "kmers")
#       snow::clusterEvalQ(cl, Rcpp::sourceCpp(path_to_rcpp_script))
#       for (i in seq_along(confirmed_genomes_paths_split)) {
#         verbose_message(paste("Working on chunk", i, "of", split), verbose)
#         output <- pbapply::pblapply(confirmed_genomes_paths_split[[i]],
#                                     convert_to_kmers, cl = cl)
#         names(output) <- confirmed_genomes_ids_split[[i]]
#         save(output, file = batch_filename(
#           output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#         remove(output)
#       }
#       snow::stopCluster(cl)
#     }
#   } else {
#     for (i in seq_along(confirmed_genomes_paths_split)) {
#       verbose_message(paste("Working on chunk", i, "of", split), verbose)
#       output <- lapply(confirmed_genomes_paths_split[[i]], convert_to_kmers,
#                        kmers = kmers,
#                        simplify = simplify,
#                        anchor= anchor,
#                        drop_n=drop_n,
#                        verbose = verbose)
#       names(output) <- confirmed_genomes_ids_split[[i]]
#       save(output, file = batch_filename(
#         output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#       remove(output)
#     }
#   }
# }
max_batch_split <- 9999 # filenames break if greater


verbose_message <- function(s, verbosity = FALSE){
  if (verbosity) message(s)
}

flat_stringset <- function(x) {
  as.character(unlist(Biostrings::readDNAStringSet(x)))
}

genome_to_kmers <- function(x,
                             kmers =3,
                             simplify =T,
                             anchor= T,
                             drop_n=T,
                             integer_index = F) {
  message(paste("Working on", x))
  x <- flat_stringset(x)
  kmers(x,
        k = kmers,
        anchor = anchor,
        simplify = simplify,
        clean_up = drop_n,
        key_as_int = integer_index)
}

genomes_to_sparse <- function(input_dir,
                              output_dir,
                              cores = 1,
                              kmers = 3,
                              n_genomes = NULL,
                              split = 1,
                              anchor = FALSE,
                              simplify = FALSE,
                              drop_n = FALSE,
                              integer_index = FALSE,
                              random_shuffle = FALSE,
                              formats = "libsvm") {
  formats <- tolower(formats)
  if ("libsvm" %in% formats & !integer_index) {
    stop("libsvm format only support integer-indexed features,
         set using integer_index")
  }

  if (split > max_batch_split) {
    stop(paste("Max number of batches =", max_batch_split))
  }

  confirmed_genome_paths <- genome_paths_from_dir(input_dir,
                                                  n_genomes = n_genomes,
                                                  ext = "fna",
                                                  random_shuffle = random_shuffle)
  confirmed_genome_ids <- strip_filename(confirmed_genome_paths)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  split_genome_paths <- split_paths(confirmed_genome_paths, split)
  split_genome_ids <- split_paths(confirmed_genome_ids, split)

  for (i in seq_along(split_genome_paths)) {
    message(paste("Working on chunk", i, "of", split))

    if (cores == 1) {
      output <- lapply(split_genome_paths[[i]], genome_to_kmers,
                       kmers = kmers,
                       simplify = simplify,
                       anchor= anchor,
                       drop_n=drop_n,
                       integer_index = integer_index)
    }

    if (cores > 1 & Sys.info()["sysname"] != "Windows") {
      output <- parallel::mclapply(split_genome_paths[[i]],
                                   genome_to_kmers,
                                   kmers = kmers,
                                   simplify = simplify,
                                   anchor= anchor,
                                   drop_n=drop_n,
                                   integer_index = integer_index,
                                   mc.cores = cores)
    }

    if (cores > 1 & Sys.info()["sysname"] == "Windows") {
      # if running on Windows, we require a few extra steps for parallel to work
      # use snow package instead
      cl <- snow::makeCluster(cores, type = 'SOCK')
      # next expose required variables and Rcpp scripts. Note that any package
      # function calls should be done without namespace (therefore e.g. Rcpp::)
      # otherwise run library(package) on all cores
      snow::clusterExport(cl, "kmers")
      for (i in seq_along(split_genome_paths)) {
        message(paste("Working on chunk", i, "of", split))
        output <- pbapply::pblapply(split_genome_paths[[i]],
                                    convert_to_kmers, cl = cl)
      }
    }
    names(output) <- split_genome_ids[[i]]


    if ("rdata" %in% formats) {
      save(output, file = batch_filename(
        output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
    }

    if ("rds" %in% formats) {
      saveRDS(output, file = batch_filename(
        output_dir, i, paste0("_k", kmers, "_kmer_data"), ".rds"))
    }

    if ("libsvm" %in% formats) {
      write_sparse_batch(
        output,
        batch_filename(
          output_dir, i, paste0("_k", kmers, "_libsvm"), ".txt"))

      meta_data <- data.frame(genome_id = names(output))
      write.csv(
        meta_data,
        batch_filename(
          output_dir, i, paste0("_k", kmers, "_meta_data"), ".csv"))
    }

    remove(output)

    if (cores > 1 & Sys.info()["sysname"] == "Windows") {
      snow::stopCluster(cl)
    }
  }
}

genomes_to_non_sparse <- function(input_dir,
                                  output_dir,
                                  cores = 1,
                                  kmers = 3,
                                  n_genomes = NULL,
                                  split = 1,
                                  anchor = FALSE,
                                  simplify = FALSE,
                                  drop_n = FALSE,
                                  integer_index = FALSE,
                                  random_shuffle = FALSE,
                                  formats = "csv") {
  formats <- tolower(formats)
  if (split > max_batch_split) {
    stop(paste("Max number of batches =", max_batch_split))
  }

  confirmed_genome_paths <- genome_paths_from_dir(input_dir,
                                                  n_genomes = n_genomes,
                                                  ext = "fna",
                                                  random_shuffle = random_shuffle)
  confirmed_genome_ids <- strip_filename(confirmed_genome_paths)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  split_genome_paths <- split_paths(confirmed_genome_paths, split)
  split_genome_ids <- split_paths(confirmed_genome_ids, split)

  for (i in seq_along(split_genome_paths)) {

    if (cores == 1) {
      output <- lapply(split_genome_paths[[i]],
                       genome_to_kmers,
                       kmers = kmers,
                       simplify = simplify,
                       anchor = anchor,
                       drop_n = drop_n,
                       integer_index = integer_index)
    }

    if (cores > 1 & Sys.info()["sysname"] != "Windows") {
      output <- parallel::mclapply(split_genome_paths[[i]],
                                   genome_to_kmers,
                                   kmers = kmers,
                                   simplify = simplify,
                                   anchor= anchor,
                                   drop_n=drop_n,
                                   integer_index = integer_index,
                                   mc.cores = cores)
    }

    if (cores > 1 & Sys.info()["sysname"] == "Windows") {
      # if running on Windows, we require a few extra steps for parallel to work
      # use snow package instead
      cl <- snow::makeCluster(cores, type = 'SOCK')
      # next expose required variables and Rcpp scripts. Note that any package
      # function calls should be done without namespace (therefore e.g. Rcpp::)
      # otherwise run library(package) on all cores
      snow::clusterExport(cl, "kmers")
      for (i in seq_along(split_genome_paths)) {
        message(paste("Working on chunk", i, "of", split))
        output <- pbapply::pblapply(split_genome_paths[[i]],
                                    convert_to_kmers, cl = cl)
      }
    }
    names(output) <- split_genome_ids[[i]]

    if ("csv" %in% formats) {
      write_non_sparse_batch(
        kmers_to_dataframe(output, k = kmers),
        batch_filename(
          output_dir,
          i,
          paste0("_k", kmers, "_kmer_data"),
          ".csv"),
        k = kmers,
        file_format = "csv",
        overwrite = TRUE)
    }

    if ("parquet" %in% formats) {
      write_non_sparse_batch(
        kmers_to_dataframe(output, k = kmers),
        batch_filename(
          output_dir,
          i,
          paste0("_k", kmers, "_kmer_data"),
          ".parquet"),
        k = kmers,
        file_format = "parquet",
        overwrite = TRUE)
    }
    remove(output)
  }
}
#
# kmer_batch <- function(input_dir,
#                        output_dir,
#                        verbose = FALSE,
#                        cores = 1,
#                        kmers = 3,
#                        n_genomes = NULL,
#                        split = 1,
#                        anchor = FALSE,
#                        simplify = FALSE,
#                        drop_n = FALSE,
#                        formats = "libsvm",
#                        integer_index = FALSE,
#                        random_shuffle = FALSE) {
#   # max_batch_split <- 9999 # filenames break if greater
#   save_formats <- tolower(unlist(strsplit(formats, split = ",")))
# #
# #   if (split > max_batch_split) {
# #     stop(paste("Max number of batches =", max_batch_split))
# #   }
# #
# #   confirmed_genomes_paths <- list.files(input_dir,
# #                                         pattern = "\\.fna$", full.names = TRUE)
# #
# #   if (random_shuffle) {
# #     verbose_message("Shuffling genomes..", verbose)
# #     confirmed_genomes_paths <- sample(confirmed_genomes_paths)
# #   }
# #
# #   confirmed_genomes_ids <- strip_filename(confirmed_genomes_paths, "fna")
# #
# #   if (length(confirmed_genomes_paths) < 1) {
# #     stop("No .fna files found in input dir")
# #   }
#
#   # n_genomes_to_process <- ifelse(is.null(n_genomes),
#   #                                length(confirmed_genomes_paths),
#   #                                n_genomes)
#   #
#   # confirmed_genomes_paths <- confirmed_genomes_paths[1:n_genomes_to_process]
#   # confirmed_genomes_ids <- confirmed_genomes_ids[1:n_genomes_to_process]
# #
# #   if (n_genomes_to_process > length(confirmed_genomes_paths)) {
# #     stop("Number of genomes to process greater than genomes available")
# #   }
#
#   # if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#   #
#   # split_indices <- parallel::splitIndices(
#   #   length(confirmed_genomes_paths), split)
#   #
#   # confirmed_genomes_paths_split <- lapply(split_indices, function(x) {
#   #   sapply(x, function(y) confirmed_genomes_paths[[y]])
#   # })
#   # confirmed_genomes_ids_split <- lapply(split_indices, function(x) {
#   #   sapply(x, function(y) confirmed_genomes_ids[[y]])
#   # })
#
#   if (cores > 1) {
#     if (Sys.info()["sysname"] != "Windows") {
#       for (i in seq_along(confirmed_genomes_paths_split))
#       {
#
#         if ("libsvm" %in% save_formats || integer_index) {
#           verbose_message(
#             paste("Working on int indexed chunk", i, "of", split),
#             verbose)
#           # libsvm format needs integer index insted of kmer string, so different
#           # call to Rcpp library
#           # also needs integer index if -i flagged
#           output <- parallel::mclapply(confirmed_genomes_paths_split[[i]],
#                                        convert_to_int_indexed_kmers,
#                                        mc.cores = cores,
#                                        simplify = simplify)
#           names(output) <- confirmed_genomes_ids_split[[i]]
#           if ("libsvm" %in% save_formats) {
#             # write_sparse_batch(
#             #   output,
#             #   batch_filename(
#             #     output_dir, i, paste0("_k", kmers, "_libsvm"), ".txt"))
#
#             # write meta data separately for libsvm
#           #   meta_data <- data.frame(genome_id = names(output))
#           #   write.csv(
#           #     meta_data,
#           #     batch_filename(
#           #       output_dir, i, paste0("_k", kmers, "_meta_data"), ".csv"))
#           # }
#
#           # if ("rdata" %in% save_formats) {
#           #   save(output, file = batch_filename(
#           #     output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#           # }
#           #
#           # if ("rds" %in% save_formats) {
#           #   rds_filepath <- batch_filename(
#           #     output_dir, i, paste0("_k", kmers, "_kmer_data"), ".rds")
#           #   verbose_message(paste("Writing rds file to", rds_filepath), verbose)
#           #   saveRDS(output, file = rds_filepath)
#           # }
#
#           remove(output)
#
#           if (integer_index) {
#             # formats below are not integer indexed
#             next()
#           }
#         }
#
#         non_integer_indexed_formats <- c("rdata", "rds", "csv", "parquet")
#         if (any(save_formats %in% non_integer_indexed_formats)) {
#           verbose_message(
#             paste("Working on chunk", i, "of", split),
#             verbose)
#
#           output <- parallel::mclapply(confirmed_genomes_paths_split[[i]],
#                                        convert_to_kmers, mc.cores = cores)
#           names(output) <- confirmed_genomes_ids_split[[i]]
#         }
#
#         if ("rdata" %in% save_formats) {
#           save(output, file = batch_filename(
#             output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#         }
#
#         if ("rds" %in% save_formats) {
#           saveRDS(output, file = batch_filename(
#             output_dir, i, paste0("_k", kmers, "_kmer_data"), ".rds"))
#         }
#
#         if ("csv" %in% save_formats || "parquet" %in% save_formats) {
#           # further processing required for these formats
#           if (!simplify) {
#             message("CSV or parquet save without --simplify is not implemented")
#           } else {
#             output <- lapply(output, function(x) {
#               x <- unlist(x)
#               x})
#
#             # keep only items with length == 4 ^ kmers
#             target_length <- 4 ^ kmers
#             output <- subset(
#               output,
#               sapply(output, function(x) {
#                 ifelse(length(x) == target_length, TRUE, FALSE)
#               }))
#
#             genome_ids <- as.character(names(output))
#             data.table::setDT(output)
#             output <- data.table::transpose(output)
#             output <- data.table(genome_id = genome_ids, output)
#             colnames(output) <- c("genome_id", paste0(
#               "kmer_", seq(from = 1, to = ncol(output) - 1)))
#             output[, genome_id := genome_ids]
#
#             if ("parquet" %in% save_formats) {
#               backup_kmer_data(
#                 output,
#                 batch_filename(
#                   output_dir,
#                   i,
#                   paste0("_k", kmers, "_kmer_data"),
#                   ".parquet"),
#                 file_format = "parquet")
#             }
#
#             if ("csv" %in% save_formats) {
#               backup_kmer_data(
#                 output,
#                 batch_filename(
#                   output_dir,
#                   i,
#                   paste0("_k", kmers, "_kmer_data"),
#                   ".csv"),
#                 overwrite = TRUE)
#             }
#
#           }
#           if (exists("output")) remove(output)
#         }
#       }
#     } else {
#       # if running on Windows, we require a few extra steps for parallel to work
#       # use snow package instead
#       cl <- snow::makeCluster(cores, type = 'SOCK')
#       # next expose required variables and Rcpp scripts. Note that any package
#       # function calls should be done without namespace (therefore e.g. Rcpp::)
#       # otherwise run library(package) on all cores
#       snow::clusterExport(cl, "kmers")
#       snow::clusterEvalQ(cl, Rcpp::sourceCpp(path_to_rcpp_script))
#       for (i in seq_along(confirmed_genomes_paths_split)) {
#         verbose_message(paste("Working on chunk", i, "of", split), verbose)
#         output <- pbapply::pblapply(confirmed_genomes_paths_split[[i]],
#                                     convert_to_kmers, cl = cl)
#         names(output) <- confirmed_genomes_ids_split[[i]]
#         save(output, file = batch_filename(
#           output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#         remove(output)
#       }
#       snow::stopCluster(cl)
#     }
#   } else {
#     for (i in seq_along(confirmed_genomes_paths_split)) {
#       verbose_message(paste("Working on chunk", i, "of", split), verbose)
#       output <- lapply(confirmed_genomes_paths_split[[i]], convert_to_kmers,
#                        kmers = kmers,
#                        simplify = simplify,
#                        anchor= anchor,
#                        drop_n=drop_n,
#                        verbose = verbose)
#       names(output) <- confirmed_genomes_ids_split[[i]]
#       save(output, file = batch_filename(
#         output_dir, i, paste0("_k", kmers, "_kmer_data"), ".RData"))
#       remove(output)
#     }
#   }
# }
