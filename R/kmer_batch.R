#' Convert genomes to kmers in libsvm format
#'
#' @param source_dir directory containing genomes
#' @param target_dir target directory to store kmers in libsvm format
#' @param k k-mer length
#' @param canonical only count canonical kmers
#' @param squeeze remove non-canonical kmers
#' @param ext file extension to filter
#'
#' @description
#' Raw genome data (pre- or post-assembly) is usually transformed by k-mer
#' counting prior to machine learning (ML). XGBoost is a popular ML algorithm
#' for this problem, due to its scalability to high dimensional data. This
#' function converts genomes to k-mer counts stored in XGBoost's preferred
#' format, libsvm. Further information on the libsvm format is available at
#' https://xgboost.readthedocs.io/en/stable/tutorials/input_format.html.
#' Briefly, libsvm is effectively a text file that stores data points as
#' x:y pairs, where x is the feature index, and y is the feature value. Each
#' observation is stored on its own line, with the first column reserved for
#' labels. Labels can be provided later, during data import.
#'
#' This function converts each individual genome to an individual libsvm
#' text file of k-mer counts (therefore, each .txt file will be 1 line long).
#' The function supports parallel processing using the future package, and
#' progress bars using the progressr package (see examples).
#'
#' Although XGBoost can load a multiple .txt (libsvm) files by providing the
#' directory as an input, this is generally not recommended as order of
#' import cannot be guaranteed and probably depends on filesystem. Instead,
#' it is recommended that this function is combined with
#' split_and_combine_files() which generates a single .txt file (with the
#' order of observations guaranteed and stored in a .csv file).
#'
#' @examples
#' \dontrun{
#'   future::plan(future::multisession)
#'   progressr::with_progress(
#'     genomes_to_kmer_libsvm("path_in", "path_out")
#'   )
#' }
#'
#' @seealso to convert a single genome, use [genome_to_libsvm()]
#' @export
genomes_to_kmer_libsvm <- function(source_dir,
                                   target_dir,
                                   k = 3,
                                   canonical = TRUE,
                                   squeeze = FALSE,
                                   ext = ".fna") {
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  ext <- gsub("^\\.", "", ext)
  genome_paths <- list.files(source_dir,
                             pattern = paste0("*.", ext),
                             full.names = TRUE,
                             ignore.case = TRUE)
  p <- progressr::progressor(along = genome_paths)
  future.apply::future_lapply(genome_paths, \(x) {
    genome_to_libsvm(as.character(Biostrings::readDNAStringSet(x)),
                    file.path(normalizePath(target_dir),
                              paste0(strip_filename(x), ".txt")),
                    k = k,
                    canonical = canonical,
                    squeeze = squeeze)
    p(glue::glue("Completed: {basename(x)}"))
  }, future.seed = TRUE)
  return(TRUE)
}
