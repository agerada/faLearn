patric_ftp_path <- "ftp://ftp.bv-brc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt"

#' Load PATRIC database
#'
#' @param x Character path to local or ftp path (.txt or .rds), or
#' data.frame object.
#'
#' @return PATRIC database (S3 class 'patric_db')
#' @export
#'
#' @examples
#' \donttest{
#' patric_db <- load_patric_db()  # will get from PATRIC ftp
#' }
#'
#' # make data.frame with single row
#' p <- data.frame(genome_id = 1,
#'                 genome_name = "E. coli",
#'                 antibiotic = "amoxicillin",
#'                 measurement = 2.0,
#'                 measurement_unit = "mg/L",
#'                 laboratory_typing_method = "Agar dilution",
#'                 resistant_phenotype = "R")
#' load_patric_db(p)
load_patric_db <- function(x = patric_ftp_path) {
  if (inherits(x, "patric_db")) {
    message("Input to load_patric_db appears to already be a patric_db")
    return(x)
  }
  if (inherits(x, "data.frame")) {
    tryCatch(check_valid_patric_db(x),
             error = function(e) stop("Data does not appear to be a compatible dataframe.
                                      See ?load_patric_db"))
    return(as_patric_db(x))
  }

  if (endsWith(x, ".rds")) {
    tryCatch(patric_db <- readRDS(x),
             error = function(e) e)
    if (inherits(patric_db, "patric_db")) {
      return(patric_db)
    } else{
      stop(".rds does not appear to be a valid patric_db. Load from .txt first
or use tidy_patric_meta_data()")
    }
  }
  if (!endsWith(x, ".txt")) {
    stop("Path to PATRIC database must be to a .txt file")
  }
  if (startsWith(x, "ftp")) {
    # Download the file using download_patric_db to handle FTPS
    temp_file <- tempfile(fileext = ".txt")
    success <- download_patric_db(temp_file, ftp_path = x)
    if (!success) {
      stop("Failed to download PATRIC database from FTP")
    }
    x <- temp_file
  }
  patric_db <- readr::read_delim(x, delim = "\t",
                                 col_types = readr::cols(.default = "c"))
  as_patric_db(patric_db)
}

check_valid_patric_db <- function(x) {
  # check inherits dataframe
  if (!inherits(x, "data.frame")) {
    stop("Data is not a data.frame")
  }

  # check columns
  required_columns <- c("genome_id", "genome_name", "antibiotic",
                        "measurement", "measurement_unit",
                        "laboratory_typing_method", "resistant_phenotype")
  if (!all(required_columns %in% colnames(x))) {
    stop("Data does not contain all required columns for PATRIC-style database.
         Please see ?load_patric_db()")
  }
}

as_patric_db <- function(x) {
  check_valid_patric_db(x)
  class(x) <- append(class(x), "patric_db", after = 0)
  x
}

#' Download PATRIC database
#'
#' @param save_path Save path (should be .txt)
#' @param ftp_path PATRIC database FTP path to download
#' @param overwrite Force overwrite
#'
#' @return TRUE if successful, FALSE if failure.

#' @export
#' @examples
#' \donttest{
#' download_patric_db(tempfile(fileext = ".txt"))
#' }
download_patric_db <- function(save_path,
                           ftp_path = patric_ftp_path,
                           overwrite = FALSE) {
  if (file.exists(save_path) & !overwrite) {
    stop("File already exists, use overwrite (carefully)")
  }
  if (!endsWith(save_path, ".txt")) {
    warning("The path provided is not a .txt path, recommend use .txt")
  }
  target_dir <- dirname(save_path)
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)

  # Try a FTPS-capable strategy. The server requires SSL/TLS on the control
  # channel (explicit FTPS / AUTH TLS). Prefer the 'curl' package (libcurl)
  # with TLS required; fall back to system curl with --ssl-reqd; finally to
  # default download.file (least reliable for FTPS).
  download_attempt <- function() {
    # Normalize legacy hostnames to the current BV-BRC FTP host (certificate matches bv-brc.org)
    ftp_path_local <- sub("ftp://ftp\\.bvbrc\\.org", "ftp://ftp.bv-brc.org", ftp_path)
    ftp_path_local <- sub("ftps://ftp\\.bvbrc\\.org", "ftps://ftp.bv-brc.org", ftp_path_local)

    if (startsWith(ftp_path_local, "ftp://") || startsWith(ftp_path_local, "ftps://")) {
      # Prefer to request the FTP URL and negotiate TLS (explicit FTPS).
      p <- ftp_path_local

      # 1) Preferred: curl package with TLS required and anonymous login.
      if (requireNamespace("curl", quietly = TRUE)) {
        handle <- curl::new_handle()
        curl::handle_setopt(handle,
                            use_ssl = 3,   # require TLS for control+data
                            userpwd = "anonymous:")
        res <- tryCatch({
          curl::curl_download(p, save_path, handle = handle, mode = "wb")
          0
        }, error = function(e) e, warning = function(w) w)
        if (is.numeric(res) && res == 0) return(0)
      }

      # 2) Fallback: system curl via utils::download.file(method = "curl")
      #    using explicit FTPS requirement and anonymous user.
      res <- tryCatch({
        utils::download.file(p, save_path, method = "curl",
                             extra = "--ssl-reqd -u anonymous:", mode = "wb")
      }, warning = function(w) w, error = function(e) e)
      if (is.numeric(res) && res == 0) return(0)

      # 3) Final fallback: try default download.file on the original path
      return(utils::download.file(ftp_path_local, save_path, mode = "wb"))
    }

    # Non-FTP(S) URLs: use the default download.file behaviour.
    return(utils::download.file(ftp_path_local, save_path, mode = "wb"))
  }

  return_val <- tryCatch(download_attempt(),
                         error = function(e) {
                           warning("Error downloading PATRIC file: ", conditionMessage(e))
                           return(1)
                         })

  if (is.numeric(return_val) && return_val != 0) {
    warning("Non-zero return value on file download")
    return(FALSE)
  }
  return(TRUE)
}

#' Automated download of genomes from PATRIC database
#'
#' @param output_directory local directory to save to
#' @param taxonomic_name character of taxonomic bacterial name to download
#' @param database local or ftp path to PATRIC database, or loaded database using load_patric_db()
#' @param filter "MIC" or "disk" or "all" phenotypes
#' @param ab antibiotic(s) of interest, provided as a character vector of
#' antibiotic names/codes, or ideally, as AMR::ab classes, created using AMR::as.ab
#' (default = all)
#' @param n_genomes number of genomes (0 = all)
#'
#' @return The number of failed downloads (i.e., 0 if all attempted downloads
#' were successful).
#'
#' @export
#' @examples
#' \donttest{
#' pull_PATRIC_genomes(tempdir(),
#'                     taxonomic_name = "Escherichia coli",
#'                     filter = "MIC",
#'                     n_genomes = 10)
#'}
pull_PATRIC_genomes <- function(output_directory,
                                taxonomic_name = NULL,
                                database = patric_ftp_path,
                                filter = "MIC",
                                ab = NULL,
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

  if (is.null(taxonomic_name)) {
    filtered_data <- patric_amr_list
  } else {
    filtered_data <- patric_amr_list[grep(taxonomic_name, patric_amr_list$genome_name), ]
  }

  filtered_data <- filtered_data |> dplyr::filter(dplyr::case_when(
    filter == "mic" & measurement_unit == "mg/L" ~ TRUE,
    filter == "disc" & laboratory_typing_method == "Disk diffusion" ~ TRUE,
    filter == "all" ~ TRUE
  ))

  if (!is.null(ab)) {
    ab <- AMR::as.ab(ab)
    filtered_data$antibiotic <- AMR::as.ab(filtered_data$antibiotic)
    filtered_data <- filtered_data |>
      dplyr::filter(.data[["antibiotic"]] %in% ab)

    # make sure mic/disk is valid
    if (filter == "mic") {
      filtered_data <- filtered_data |>
        dplyr::mutate(measurement = AMR::as.mic(MIC::clean_raw_mic(.data[["measurement"]]))) |>
        dplyr::filter(!is.na(.data[["measurement"]]))
    } else if (filter == "disc") {
      filtered_data <- filtered_data |>
        dplyr::mutate(measurement = AMR::as.disk(.data[["measurement"]])) |>
        dplyr::filter(!is.na(.data[["measurement"]]))
    }

    # if filter == "all", then all measurements are kept (even if turn out to
    # to be invalid).
  }

  genome_ids <- unique(filtered_data$genome_id)

  if (n_genomes < 0) {
    n_downloads <- 0
  } else if (n_genomes > 0 & n_genomes < length(genome_ids)) {
    n_downloads <- n_genomes
  } else {
    n_downloads <- length(genome_ids)
  }

  genome_paths <- glue::glue(
    "ftp://ftp.bv-brc.org/genomes/{genome_ids}/{genome_ids}.fna"
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
      # Normalize legacy hostnames
      p <- genome_paths[[i]]
      p_local <- sub("ftp://ftp\\.patricbrc\\.org", "ftp://ftp.bv-brc.org", p)
      p_local <- sub("ftps://ftp\\.patricbrc\\.org", "ftps://ftp.bv-brc.org", p_local)

      downloaded <- FALSE

      # 1) Preferred: curl package with TLS required and anonymous login.
      if (requireNamespace("curl", quietly = TRUE)) {
        handle <- curl::new_handle()
        curl::handle_setopt(handle, use_ssl = 3, userpwd = "anonymous:")
        res <- tryCatch({
          curl::curl_download(p_local, target_path, handle = handle, mode = "wb")
          TRUE
        }, error = function(e) FALSE, warning = function(w) FALSE)
        downloaded <- isTRUE(res)
      }

      # 2) Fallback: system curl via utils::download.file(method = "curl") with FTPS
      if (!downloaded) {
        res2 <- tryCatch({
          utils::download.file(p_local, destfile = target_path, method = "curl",
                               extra = "--ssl-reqd -u anonymous:", mode = "wb")
        }, warning = function(w) w, error = function(e) e)
        downloaded <- is.numeric(res2) && res2 == 0
      }

      # 3) Final fallback: default download.file
      if (!downloaded) {
        res3 <- tryCatch({
          utils::download.file(p_local, destfile = target_path, mode = "wb")
        }, warning = function(w) w, error = function(e) e)
        downloaded <- is.numeric(res3) && res3 == 0
      }

      if (!downloaded) {
        failures <- failures + 1
        message(glue::glue("Unable to download {genome_ids[[i]]}"))
      }
    }
    i <- i + 1
  }
  failures
}

#' Tidy PATRIC data
#'
#' @param x PATRIC database loaded using MIC::load_patric_db
#' @param prefer_more_resistant High MICs, narrow zones, or resistant phenotypes
#' will be preferred where multiple reported for the same isolate
#' @param as_ab convert antibiotics to AMR::ab class (column names are antibiotic
#' codes)
#' @param filter_abx filter antibiotics of interest, provided as a vector of
#' antibiotics character names/codes, or ideally, as AMR::ab classes, created
#' using AMR::as.ab
#'
#' @return Tidy data, with antimicrobials in wide format, column names describing
#' methodology ("mic_", "disk_", "pheno_"). S3 class "tidy_patric_db".
#' @export
#' @examples
#' db <- data.frame(genome_id = 1,
#'                 genome_name = "E. coli",
#'                 antibiotic = "amoxicillin",
#'                 measurement = 2.0,
#'                 measurement_unit = "mg/L",
#'                 laboratory_typing_method = "Agar dilution",
#'                 resistant_phenotype = "R")
#' db <- load_patric_db(db)
#' tidy_patric_meta_data(db)
tidy_patric_meta_data <- function(x,
                                  prefer_more_resistant = TRUE,
                                  as_ab = TRUE,
                                  filter_abx = NULL) {
  if (!inherits(x, "patric_db")) {
    stop("Please load data using MIC::load_patric_db()")
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
    dplyr::filter(.data[["laboratory_typing_method"]] %in% c("Agar dilution",
                                                             "Broth dilution",
                                                             "MIC")) |>
    dplyr::mutate(measurement = AMR::as.mic(MIC::clean_raw_mic(.data[["measurement"]]))) |>
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

  output <- dplyr::full_join(output, pheno_data, by = c("genome_id", "genome_name"))
  output <- as.data.frame(output)
  class(output) <- append(class(output), "patric_db", after = 0)
  class(output) <- append(class(output), "tidy_patric_db", after = 0)
  return(output)
}
