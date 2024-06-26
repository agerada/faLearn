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
  suppressWarnings(
    dplyr::case_when(
      stringr::str_detect(mic, ">") ~ AMR::as.mic(mic * scale),
      stringr::str_detect(mic, "<") ~ AMR::as.mic(mic / scale),
      .default = mic)
  )
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
    return(utils::head(start, -1))
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
#' @param levels_from_AMR conform to AMR::as.mic levels
#' @param max_conc maximum concentration to force to
#' @param min_conc minimum concentration to force to
#' @param prefer where value is in between MICs (e.g., 24mg/L) chose the higher
#' MIC ("max") or lower MIC ("min").
#'
#' @return AMR::as.mic compatible character
#' @export
#'
force_mic <- function(value,
                      levels_from_AMR = TRUE,
                      max_conc = 2048,
                      min_conc = 0.00005,
                      prefer = 'max') {

  if (levels_from_AMR) {
    mic_levels <- levels(AMR::as.mic(NA))
  } else {
    mic_levels <- mic_range(start = max_conc, min = min_conc)
  }

  output <- rep(NA_character_, length(value))
  for (i in seq_along(value)) {

    prefix <- NULL
    inner_val <- value[i]
    if (AMR::is.mic(inner_val)) {
      inner_val <- as.character(inner_val)
    }
    if (is.numeric(inner_val)) {
      appropriate_levels <- subset(mic_levels, !stringr::str_detect(mic_levels, "[^0-9.]"))
    }
    if (is.character(inner_val)){
      if (stringr::str_detect(inner_val, "<")) {
        appropriate_levels <- subset(mic_levels, stringr::str_detect(mic_levels, "<"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        prefix <- "<"
      } else if (stringr::str_detect(inner_val, ">")) {
        appropriate_levels <- subset(mic_levels, stringr::str_detect(mic_levels, ">"))
        appropriate_levels <- stringr::str_remove_all(appropriate_levels, "[^0-9.]")
        prefix <- ">"
      } else {
        appropriate_levels <- subset(mic_levels, !stringr::str_detect(mic_levels, "[^0-9.]"))
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

essential_agreement <- function(gold_standard, test) {
  frac <- test / gold_standard
  return(factor(
    dplyr::case_when(
      frac > 2.0 ~ FALSE,
      frac < 0.5 ~ FALSE,
      frac != 1.0 ~ TRUE,
      TRUE ~ TRUE)))
}

compare_mic <- function(gold_standard, test, ab = NA, simplify = TRUE) {
  gold_standard <- force_mic(gold_standard)
  test <- force_mic(test)

  gold_standard <- as.character(gold_standard)
  test <- as.character(test)

  gold_standard <- mic_uncensor(gold_standard)
  test <- mic_uncensor(test)

  gold_standard <- as.numeric(gold_standard)
  test <- as.numeric(test)

  output <- data.frame(
    ab = ab,
    gold_standard = AMR::as.mic(gold_standard),
    test = AMR::as.mic(test),
    essential_agreement = essential_agreement(gold_standard, test)
  )
  class(output) <- append(class(output), "mic_validation", 0)
  output
}

compare_sir <- function(gold_standard, test) {
  return(
    factor(
      dplyr::case_when(
        gold_standard == "S" & test == "I" ~ "minor",
        gold_standard == "R" & test == "I" ~ "minor",
        gold_standard == "I" & test == "S" ~ "minor",
        gold_standard == "I" & test == "R" ~ "minor",
        gold_standard == "S" & test == "R" ~ "major",
        gold_standard == "R" & test == "S" ~ "very major",
        TRUE ~ NA
      )
    )
  )
}

#' Plot MIC validation results
#'
#' @param x object generated using compare_mic
#' @param ... additional arguments
#'
#' @export
plot.mic_validation <- function(x, ...) {
  x |>
    dplyr::group_by(.data[["gold_standard"]],
                    .data[["test"]],
                    .data[["essential_agreement"]],
                    .data[["ab"]]) |>
    dplyr::mutate(ab = AMR::ab_name(AMR::as.ab(as.character(.data[["ab"]])))) |>
    dplyr::mutate(ab = dplyr::if_else(is.na(.data[["ab"]]), "unknown", .data[["ab"]])) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::rename(`EA` = .data[["essential_agreement"]]) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[["gold_standard"]],
                                 y = .data[["test"]],
                                 fill = .data[["n"]],
                                 color = .data[["EA"]])) +
    ggplot2::geom_tile(alpha=1) +
    ggplot2::geom_text(ggplot2::aes(label=.data[["n"]])) +
    ggplot2::scale_fill_gradient(low="white", high="#009194") +
    ggplot2::scale_fill_manual(values=c("red", "black"), aesthetics = "color") +
    lemon::facet_rep_wrap(~ .data[["ab"]], nrow = 2, repeat.tick.labels = TRUE) +
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(fill=NA))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::theme(legend.position = c(0.9,0.2)) +
    ggplot2::xlab("Gold standard MIC (mg/L)") +
    ggplot2::ylab("Test (mg/L)")
}
