#' Get MIC meta-data from feature database
#'
#' @param x dataframe containing meta-data
#' @param ids vector of IDs to get meta-data for
#' @param id_col column name containing IDs
#' @param ab_col column name containing MIC results
#' @param as_mic return as AMR::as.mic
#' @param prefer_high_mic where multiple MIC results per ID, prefer the higher MIC
#' @param simplify return as vector of MICs (vs dataframe)
#'
#' @return vector containing MICs, or dataframe of IDs and MICs
#' @export
#'
#' @examples
#' df <- data.frame(genome_id = c("a_12", "b_42", "x_21", "x_21", "r_75"),
#'                  gentamicin = c(0.25, 0.125, 32.0, 16.0, "<0.0125"))
#' get_mic(df,
#'         ids = c("b_42", "x_21"),
#'         id_col = "genome_id",
#'         ab_col = "gentamicin",
#'         as_mic = FALSE,
#'         prefer_high_mic = TRUE,
#'         simplify = TRUE)
get_mic <- function(x,
                    ids,
                    id_col,
                    ab_col,
                    as_mic = TRUE,
                    prefer_high_mic = TRUE,
                    simplify = TRUE) {
  if ("save_order" %in% names(x)) {
    stop("Unable to work with x that contains column name 'save_order', please
         rename.")
  }

  x <- x[order(x[[id_col]], x[[ab_col]], decreasing = prefer_high_mic),]
  x <- x[!duplicated(x[[id_col]]),]

  output <- data.frame(ids)
  names(output) <- id_col
  output$save_order <- 1:nrow(output)
  output <- merge(output, x,
                  by = id_col,
        sort = FALSE,
        all.x = TRUE)
  if (as_mic) {
    rlang::check_installed("AMR", "To return as MIC class, AMR package must be installed")
    output[[ab_col]] <- AMR::as.mic(output[[ab_col]])
  }
  output <- output[order(output$save_order),]
  output$save_order <- NULL
  if (simplify) {
    return(output[[ab_col]])
  } else {
    return(output)
  }
}

#' Clean up raw MIC for use as a feature
#'
#' @param mic character containing MIC/s
#'
#' @return character of clean MIC/s
#' @description
#' Removes leading "=" which can sometimes be present in raw MIC results. Also converts co-trimoxazole to trimethprim component only.
#'
#' @export
#'
#' @examples
#' clean_raw_mic(c("==>64","0.25/8.0"))
clean_raw_mic <- function(mic) {
  mic <- stringr::str_remove(mic, pattern = "^=+")
  stringr::str_remove(mic, "/(.)+")
}

#' Uncensor MICs
#'
#' @param mic vector of MICs to uncensor; will be coerced to MIC using AMR::as.mic
#' @param scale scalar to multiply or divide MIC by
#'
#' @return vector of MICs in AMR::mic format
#' @description
#' Censored MIC data is generally unsuitable for modelling without some
#' conversion of censored data. This function halves MICs under the limit of
#' detection (<=) and doubles MICs above the limit.
#' @export
#'
#' @examples
#' mic_uncensor(c(">64.0", "<0.25", "8.0"))
mic_uncensor <- function(mic, scale = 2) {
  rlang::check_installed("AMR", "mic_uncensor requires AMR package")
  if (!inherits(mic, "mic")) {
    mic <- AMR::as.mic(mic)
  }
  dplyr::case_when(
    stringr::str_detect(mic, ">") ~ AMR::as.mic(mic * scale),
    stringr::str_detect(mic, "<") ~ AMR::as.mic(mic / scale),
    .default = mic)
}

#' Generate dilution series
#'
#' @param start starting (highest) concentration
#' @param dilutions number of dilutions
#' @param min minimum (lowest) concentration
#'
#' @return Vector of numeric concentrations
#' @export
#'
#' @examples
#' mic_range(128)
#' mic_range(128, dilutions = 21) # same results
mic_range <- function(start, dilutions = Inf, min = 0.0001) {
  if (start[length(start)] < min) {
    return(head(start, -1))
  }
  if (dilutions == 0) {
    return (start)
  } else {
    return(mic_range(c(start, start[length(start)] / 2),
                     dilutions - 1,
                     min))
  }
}

#' Force MIC-like into MIC-compatible format
#'
#' @param value vector of MIC-like values (numeric or character)
#' @param max_conc
#' @param min_conc
#' @param prefer
#'
#' @return
#' @export
#'
#' @examples
force_mic_again <- function(value, max_conc = 2048, min_conc = 0.00005, prefer = 'max') {
  amr_levels <- levels(AMR::as.mic(NA))
  output <- rep(NA_character_, length(value))
  for (i in seq_along(value)) {

    prefix <- NULL
    inner_val <- value[i]
    if (is.numeric(inner_val)) {
      appropriate_levels <- subset(amr_levels, !stringr::str_detect(amr_levels, "[^0-9.]"))
    }
    if (is.character(inner_val)){
      if (stringr::str_detect(inner_val, "<")) {
        appropriate_levels <- subset(amr_levels, stringr::str_detect(amr_levels, "<"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        prefix <- "<"
      } else if (stringr::str_detect(inner_val, ">")) {
        appropriate_levels <- subset(amr_levels, stringr::str_detect(amr_levels, ">"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        prefix <- ">"
      } else {
        appropriate_levels <- subset(amr_levels, !stringr::str_detect(amr_levels, "[^0-9.]"))
      }
      inner_val <- stringr::str_remove_all(inner_val, "[^0-9.]")
      inner_val <- as.numeric(inner_val)
    }

    mic_vector <- sort(as.numeric(appropriate_levels))
    positions <- which(abs(mic_vector - inner_val) == min(abs(mic_vector - inner_val)))
    if (length(positions) == 1) {
      output[i] <- paste0(prefix, mic_vector[positions])
    }
    if (prefer == 'min') {
      output[i] <- paste0(prefix, mic_vector[min(positions)])
    }
    output[i] <- paste0(prefix, mic_vector[max(positions)])
  }
  return(output)
}
