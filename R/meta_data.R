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
  pre_sort_ids <- data.frame(id = x[[id_col]])
  x <- x[order(x[[id_col]], x[[ab_col]], decreasing = prefer_high_mic),]
  x <- x[!duplicated(x[[id_col]]),]

  output <- data.frame(ids)
  names(output) <- id_col
  output <- merge(output, x,
        sort = FALSE,
        all.x = TRUE)
  if (as_mic) {
    rlang::check_installed("AMR", "To return as MIC class, AMR package must be installed")
    output[[ab_col]] <- AMR::as.mic(output[[ab_col]])
  }
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
