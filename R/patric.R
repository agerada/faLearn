patric_ftp_path <- "ftp://ftp.bvbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt"

#' Load PATRIC database
#'
#' @param path Character to local or ftp path (.txt file)
#'
#' @return PATRIC database (S3 class 'patric_db')
#' @export
#'
#' @examples
#' \dontrun{
#' patric_db <- load_patric_db()
#' }
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

#' Save PATRIC database locally
#'
#' @param save_path Save path (should be .txt)
#' @param ftp_path PATRIC database FTP path to download
#' @param overwrite Force overwrite
#'
#' @export
#'
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

#' Automated download of genomes from PATRIC database
#'
#' @param database local or ftp path to PATRIC database, or loaded database using load__patric_db()
#' @param taxonomic_name character of taxonomic bacterial name to download
#' @param filter "MIC" or "disk" or "all" phenotypes
#' @param output_directory local directory to save to
#' @param n_genomes number of genomes (0 = all)
#'
#' @export
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

#' Tidy PATRIC data
#'
#' @param x PATRIC database loaded using molMIC::load_patric_db
#' @param prefer_more_resistant High MICs, narrow zones, or resistant phenotypes
#' will be preferred where multiple reported for the same isolate
#' @param as_ab convert antibiotics to AMR::ab class (column names are antibiotic
#' codes)
#' @param filter_abx filter antibiotics of interest, provided as a vector of
#' antibiotics character names/codes, or ideally, as AMR::ab classes, created
#' using AMR::as.ab
#'
#' @return Tidy data, with antimicrobials in wide format, column names describing
#' methodology ("mic_", "disk_", "pheno_")
#' @export
tidy_patric_meta_data <- function(x,
                                  prefer_more_resistant = TRUE,
                                  as_ab = TRUE,
                                  filter_abx = NULL) {
  if (!inherits(x, "patric_db")) {
    stop("Please load data using molMIC::load_patric_db()")
  }

  if (isTRUE(as_ab) | !is.null(filter_abx)) {
    rlang::check_installed("AMR", "Antibiotic-specific arguments need AMR package")
  }

  if (as_ab) {
    x$antibiotic <- AMR::as.ab(x$antibiotic)
  }

  if (!is.null(filter_abx)) {
    filter_abx <- AMR::as.ab(filter_abx)
    x <- dplyr::filter(x, .data[["antibiotic"]] %in% filter_abx)
  }

  aggregate_mic <- list(which.min, which.max)
  mic_data <- x |>
    dplyr::filter(.data[["laboratory_typing_method"]] %in% c("Agar dilution", "Broth dilution")) |>
    dplyr::mutate(measurement = AMR::as.mic(clean_raw_mic(.data[["measurement"]]))) |>
    dplyr::group_by(.data[["genome_id"]], .data[["antibiotic"]]) |>
    dplyr::slice(aggregate_mic[[prefer_more_resistant + 1]](.data[["measurement"]])) |>
    tidyr::pivot_wider(id_cols = c(.data[["genome_id"]], .data[["genome_name"]]),
                names_from = .data[["antibiotic"]],
                values_from = .data[["measurement"]],
                names_prefix = "mic_")

  aggregate_disk <- list(which.max, which.min)
  disk_data <- x |>
    dplyr::filter(.data[["laboratory_typing_method"]] %in% c("Disk diffusion")) |>
    dplyr::mutate(measurement = AMR::as.disk(.data[["measurement"]])) |>
    dplyr::group_by(.data[["genome_id"]], .data[["antibiotic"]]) |>
    dplyr::slice(aggregate_disk[[prefer_more_resistant + 1]](.data[["measurement"]])) |>
    tidyr::pivot_wider(id_cols = c(.data[["genome_id"]], .data[["genome_name"]]),
                       names_from = .data[["antibiotic"]],
                       values_from = .data[["measurement"]],
                       names_prefix = "disk_")
  output <- dplyr::full_join(mic_data, disk_data, by = c("genome_id", "genome_name"))

  aggregate_sir <- list(which.min, which.max)
  pheno_data <- x |>
    dplyr::filter(!is.na(.data[["resistant_phenotype"]])) |>
    dplyr::mutate(resistant_phenotype = AMR::as.sir(.data[["resistant_phenotype"]])) |>
    dplyr::group_by(.data[["genome_id"]], .data[["antibiotic"]]) |>
    dplyr::slice(aggregate_sir[[prefer_more_resistant + 1]](.data[["resistant_phenotype"]])) |>
    tidyr::pivot_wider(id_cols = c(.data[["genome_id"]], .data[["genome_name"]]),
                       names_from = .data[["antibiotic"]],
                       values_from = .data[["resistant_phenotype"]],
                       names_prefix = "pheno_")

  dplyr::full_join(output, pheno_data, by = c("genome_id", "genome_name"))
}
