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
                             integer_index = F,
                            starting_index = 1) {
  message(paste("Working on", x))
  x <- flat_stringset(x)
  kmers(x,
        k = kmers,
        anchor = anchor,
        simplify = simplify,
        clean_up = drop_n,
        key_as_int = integer_index,
        starting_index = starting_index)
}

#' Convert genomes to kmer dataset
#'
#' @param input_dir directory containing .fna files
#' @param output_dir directory to store output
#' @param cores number of cores for parallel processing (1 for sequential)
#' @param kmers kmer length
#' @param n_genomes number of genomes to process (default all)
#' @param split split disk storage of dataset into multiple files
#' @param anchor kmer algorithm parameter (see ?molMIC::kmers)
#' @param simplify kmer algorithm parameter (see ?molMIC::kmers)
#' @param drop_n kmer algorithm parameter (see ?molMIC::kmers)
#' @param integer_index kmer algorithm parameter (see ?molMIC::kmers)
#' @param starting_index kmer algorithm parameter (see ?molMIC::kmers)
#' @param random_shuffle randomise dataset
#' @param formats character vector of formats. Supports libsvm (default), parquet, csv, rds, rdata.
#'
#' @export
genomes_to_kmer_dataset <- function(input_dir,
                              output_dir,
                              cores = 1,
                              kmers = 3,
                              n_genomes = NULL,
                              split = 1,
                              anchor = FALSE,
                              simplify = FALSE,
                              drop_n = TRUE,
                              integer_index = TRUE,
                              starting_index = 1,
                              random_shuffle = TRUE,
                              formats = "libsvm") {
  formats <- tolower(formats)
  if ("libsvm" %in% formats & !integer_index) {
    stop("libsvm format only support integer-indexed features,
         set using integer_index")
  }
  if ("parquet" %in% formats & !anchor) {
    warning("Parquet format is dense, only compatible with anchor = TRUE")
    return()
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
                       integer_index = integer_index,
                       starting_index = starting_index)
    }

    if (cores > 1 & Sys.info()["sysname"] != "Windows") {
      output <- parallel::mclapply(split_genome_paths[[i]],
                                   genome_to_kmers,
                                   kmers = kmers,
                                   simplify = simplify,
                                   anchor= anchor,
                                   drop_n=drop_n,
                                   integer_index = integer_index,
                                   starting_index = starting_index,
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
                                    genome_to_kmers, cl = cl,
                                    kmers = kmers,
                                    simplify = simplify,
                                    anchor= anchor,
                                    drop_n=drop_n,
                                    integer_index = integer_index,
                                    starting_index = starting_index)
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
      utils::write.csv(
        meta_data,
        batch_filename(
          output_dir, i, paste0("_k", kmers, "_meta_data"), ".csv"))
    }

    if ("csv" %in% formats) {
      if (simplify) {
        warning("csv does not support simplify, must be a non-sparse matrix")
      } else {
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
    }

    if ("parquet" %in% formats) {
      if (simplify) {
        warning("parquet does not support simplify, must be a non-sparse matrix")
      } else {
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
    }

    remove(output)

    if (cores > 1 & Sys.info()["sysname"] == "Windows") {
      snow::stopCluster(cl)
    }
  }
}

genomes_to_kmer_libsvm <- function(source_dir,
                                   target_dir,
                                   k = 3,
                                   ext = ".fna",
                                   cores = 1) {
  ext <- gsub("^\\.", "", ext)
  genome_paths <- list.files(source_dir,
                             pattern = paste0("*.", ext),
                             full.names = TRUE,
                             ignore.case = TRUE)
  if (cores > 1) {
    parallel::mclapply(genome_paths, \(x) {
      kmers_to_libsvm(flat_stringset(x),
                      file.path(target_dir,
                                paste0(strip_filename(x), ".txt")),
                      k = k)
    })
    return(NULL)
  }
  lapply(genome_paths, \(x) {
    kmers_to_libsvm(flat_stringset(x),
                   file.path(target_dir,
                             paste0(strip_filename(x), ".txt")),
                   k = k)
  })
}
