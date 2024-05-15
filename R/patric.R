patric_ftp_path <- "ftp://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt"

#' Load PATRIC database
#'
#' @param path Character to local or ftp path (.txt file)
#'
#' @return PATRIC database (S3 class 'patric_db')
#' @export
#'
#' @examples
#' patric_db <- load_patric_db()
load_patric_db <- function(
    path = patric_ftp_path) {
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

save_patric_db <- function(save_path,
                           ftp_path = patric_ftp_path,
                           overwrite = FALSE) {
  if (file.exists(save_path) & !overwrite) {
    stop("File already exists, use overwrite (carefully)")
  }
  if (!endsWith(save_path, ".txt")) {
    warning("The path provided is not a .txt path, recommend use .txt")
  }
  target_dir <- dirname(save_path)
  if (!dir.exists(target_dir)) dir.create(target_dir)

  return_val <- utils::download.file(ftp_path, save_path, mode = "wb")
  if (return_val != 0) warning("Non-zero return value on file download")
}

pull_PATRIC_genomes <- function(database = patric_ftp_path,
                                taxonomic_name,
                                filter = "MIC",
                                output_directory,
                                n_genomes = 0) {
  supported_modality_filters <- c("all", "mic", "disc")
  filter <- tolower(filter)
  filter <- ifelse(filter == "disk", "disc", filter)

  if (!filter %in% supported_modality_filters) {
    stop(glue::glue("Unable to recognise filter {filter}, please use one of:
    {glue_collapse(supported_modality_filters, sep=', ')}"))
  }

  if (inherits(database, "patric_db")) {
    patric_amr_list <- database
  } else {
    patric_amr_list <- load_patric_db(database)
  }

  filtered_data <- patric_amr_list[grep(taxonomic_name, patric_amr_list$genome_name), ]

  filtered_data <- filtered_data |> dplyr::filter(dplyr::case_when(
    filter == "mic" & measurement_unit == "mg/L" ~ TRUE,
    filter == "disc" & laboratory_typing_method == "Disk diffusion" ~ TRUE,
    filter == "all" ~ TRUE
  ))

  genome_ids <- unique(filtered_data$genome_id)

  if (n_genomes < 0) {
    n_downloads <- 0
  } else if (n_genomes > 0 & n_genomes < length(genome_ids)) {
    n_downloads <- n_genomes
  } else {
    n_downloads <- length(genome_ids)
  }

  genome_paths <- glue::glue(
    "ftp://ftp.patricbrc.org/genomes/{genome_ids}/{genome_ids}.fna"
  )

  if (!dir.exists(output_directory)) dir.create(output_directory)

  i <- 1
  failures <- 0
  while (i <= n_downloads) {
    target_path <- file.path(output_directory,
                             glue::glue("{genome_ids[[i]]}.fna"))
    if (file.exists(target_path)) {
      message(glue::glue("Genome {genome_paths[[i]]} already exists"))
    } else {
      message(glue::glue("Downloading file {i} of {n_downloads}"))
      tryCatch(utils::download.file(genome_paths[[i]],
                             destfile = target_path,
                             mode = "wb"),
               error = function(e) {
                 failures <- failures + 1
                 message(glue::glue("Unable to download {genome_ids[[i]]}"))
               }
      )
    }
    i <- i + 1
  }
}
