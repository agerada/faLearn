load_patric_db <- function(
    path = "ftp://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt") {
  if (!endsWith(path, ".txt")) {
    stop("Path to PATRIC database must be to a .txt file")
  }
  if (startsWith(path, "ftp")) {
    path <- url(path)
  }
  patric_db <- readr::read_delim(path, delim = "\t",
                                 col_types = readr::cols(.default = "c"))
  class(patric_db) <- append(class(patric_db), "patric_db", after = 0)
  patric_db
}

pull_PATRIC_genomes <- function(database,
                                taxonomic_name,
                                filter,
                                output_directory,
                                n_genomes) {
  supported_modality_filters <- c("all", "mic", "disc")
  filter <- tolower(filter)
  filter <- ifelse(filter == "disk", "disc", filter)

  if (!filter %in% supported_modality_filters) {
    stop(glue::glue("Unable to recognise filter {filter}, please use one of:
    {glue_collapse(supported_modality_filters, sep=', ')}"))
  }
}
