
#' Convert genomes to kmers (libsvm format)
#'
#' @param source_dir directory containing genomes
#' @param target_dir target directory to store kmers in libsvm format
#' @param k kmer count
#' @param canonical only count canonical kmers
#' @param squeeze remove non-canonical kmers
#' @param ext file extension to filter
#'
#' @description
#' Supports progressr (with_progress) and future.
#'
#' @examples
#' \dontrun{
#'   future::plan(future::multisession)
#'   progressr::with_progress(
#'     genomes_to_kmer_libsvm("path_in", "path_out")
#'   )
#' }
#'
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
    kmers_to_libsvm(as.character(Biostrings::readDNAStringSet(x)),
                    file.path(normalizePath(target_dir),
                              paste0(strip_filename(x), ".txt")),
                    k = k,
                    canonical = canonical,
                    squeeze = squeeze)
    p(glue::glue("Completed: {basename(x)}"))
  }, future.seed = TRUE)
  return(TRUE)
}

#' Reverse complement of string
#'
#' @param x string
#'
#' @return reverse complement of string
#' @export
rev_comp <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
}
