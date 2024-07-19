#' Get MIC meta-data from feature database
#'
#' @param x dataframe containing meta-data
#' @param ids vector of IDs to get meta-data for
#' @param ab_col column name containing MIC results
#' @param id_col column name containing IDs
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
#'         ab_col = "gentamicin",
#'         id_col = "genome_id",
#'         as_mic = FALSE,
#'         prefer_high_mic = TRUE,
#'         simplify = TRUE)
get_mic <- function(x,
                    ids,
                    ab_col,
                    id_col = NULL,
                    as_mic = TRUE,
                    prefer_high_mic = TRUE,
                    simplify = TRUE) {
  if ("save_order" %in% names(x)) {
    stop("Unable to work with x that contains column name 'save_order', please
         rename.")
  }

  if (inherits(x, "tidy_patric_db")) {
    id_col = "genome_id"
  }

  if (!inherits(x, "tidy_patric_db") & is.null(id_col)) {
    stop("Provide id_col or pre-process meta data using tidy_patric_db()")
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

#' Calculate sample weights
#'
#' @param x Vector of discrete variables
#'
#' @return sample weights
#' @export
sample_weights <- function(x) {
  freq_table <- table(x)
  as.vector(1 / freq_table[as.character(x)])
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

  suppressWarnings(
    dplyr::case_when(
      stringr::str_detect(mic, ">") ~ AMR::as.mic(AMR::as.mic(stringr::str_remove(force_mic(mic), ">")) * scale),
      stringr::str_detect(mic, "<") ~ AMR::as.mic(AMR::as.mic(stringr::str_remove(force_mic(mic), "<")) / scale),
      .default = AMR::as.mic(force_mic(mic)))
  )
}

censor_rules <- list("B_ESCHR_COLI" = list(
  "AMK" = list(min = 2, max = 32),
  "CHL" = list(min = 4, max = 64),
  "GEN" = list(min = 1, max = 16),
  "CIP" = list(min = 0.015, max = 4),
  "MEM" = list(min = 0.016, max = 16),
  "AMX" = list(min = 2, max = 64),
  "AMC" = list(min = 2, max = 64),
  "FEP" = list(min = 0.5, max = 64),
  "CAZ" = list(min = 1, max = 128),
  "TGC" = list(min = 0.25, max = 1)
))

#' Censor MIC values
#'
#' @param mic MIC measurements
#' @param ab antibiotic name
#' @param mo microorganism name
#' @param rules censor rules
#'
#' @return cencored MIC values
#' @export
mic_censor <- Vectorize(
    function(mic, ab, mo, rules = NULL) {
      if (is.null(rules)) {
        message("No censor rules provided, using default rules")
        rules <- censor_rules
      }
      mic <- AMR::as.mic(mic)
      min_thresh <- censor_rules[[mo]][[ab]][["min"]]
      min_thresh <- ifelse(is.null(min_thresh),
                           -Inf,
                           min_thresh)
      max_thresh <- censor_rules[[mo]][[ab]][["max"]]
      max_thresh <- ifelse(is.null(max_thresh),
                           Inf,
                           max_thresh)
      if (mic > max_thresh) {
        return(paste0(">", max_thresh))
      }
      if (mic < min_thresh) {
        return(paste0("<=", min_thresh))
      }
      return(mic)
  },
  USE.NAMES = FALSE
)

#' Generate dilution series
#'
#' @param start starting (highest) concentration
#' @param dilutions number of dilutions
#' @param min minimum (lowest) concentration
#' @param precise force range to be high precision (not usually desired behaviour)
#'
#' @return Vector of numeric concentrations
#' @export
#'
#' @examples
#' mic_range(128)
#' mic_range(128, dilutions = 21) # same results
mic_range <- function(start = 512,
                      dilutions = Inf,
                      min = 0.002,
                      precise = FALSE) {
  recursive_inner <- function(start,
                              dilutions,
                              min) {
    if (start[length(start)] < min) {
      return(utils::head(start, -1))
    }
    if (dilutions == 0) {
      return (start)
    } else {
      return(recursive_inner(c(start, start[length(start)] / 2),
                       dilutions - 1,
                       min))
    }
  }

  if (precise) {
    return(recursive_inner(start = start,
                           dilutions = dilutions,
                           min = min))
  }

  eucast_range <- c(512,
                    256,
                    128,
                    64,
                    32,
                    16,
                    8,
                    4,
                    2,
                    1,
                    0.5,
                    0.25,
                    0.125,
                    0.06,
                    0.03,
                    0.016,
                    0.008,
                    0.004,
                    0.002)
  filtered_range <- eucast_range[eucast_range >= min & eucast_range <= start]
  end_index <- ifelse(is.infinite(dilutions), length(filtered_range), dilutions)
  return(filtered_range[1:end_index])
}

#' Force MIC-like into MIC-compatible format
#'
#' @param value vector of MIC-like values (numeric or character)
#' @param levels_from_AMR conform to AMR::as.mic levels
#' @param max_conc maximum concentration to force to
#' @param min_conc minimum concentration to force to
#' @param prefer where value is in between MIC (e.g., 24mg/L) chose the higher
#' MIC ("max") or lower MIC ("min").
#'
#' @return AMR::as.mic compatible character
#' @export
#'
#' @examples
#' force_mic(c("2.32", "<4.12", ">1.01"))
force_mic <- function(value,
                      levels_from_AMR = FALSE,
                      max_conc = 512,
                      min_conc = 0.002,
                      prefer = 'max') {

  if (levels_from_AMR) {
    mic_levels <- levels(AMR::as.mic(NA))
  } else {
    mic_levels <- mic_range(start = max_conc, min = min_conc)
    mic_levels <- c(max(mic_levels) * 2,
               mic_levels,
               min(mic_levels) / 2)
    mic_levels <- c(mic_levels,
                    paste0(">", mic_levels),
                    paste0("<", mic_levels))
  }

  output <- rep(NA_character_, length(value))
  for (i in seq_along(value)) {

    prefix <- NULL
    inner_val <- value[i]

    if (is.na(inner_val)) {
      next
    }

    stopifnot(any(c(AMR::is.mic(inner_val), is.numeric(inner_val), is.character(inner_val), is.na(inner_val))))

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

    if (inner_val > max_conc * 2) {
      warning(paste("Value greater than supported MIC:", inner_val))
      output[i] <- NA_character_
    } else if (inner_val < min_conc / 2) {
      print('here')
      warning(paste("Value less than supported MIC:", inner_val))
      output[i] <- NA_character_
    } else {
      output[i] <- paste0(prefix, mic_vector[max(positions)])
    }
  }
  return(output)
}

#' Essential agreement for MIC validation
#'
#' @param x AMR::mic or coercible
#' @param y AMR::mic or coercible
#' @param coerce_mic convert to AMR::mic
#' @param mode Categorical or numeric
#' @return logical vector
#' @export
#'
#' @examples
#' x <- AMR::as.mic(c("<0.25", "8", "64", ">64"))
#' y <- AMR::as.mic(c("<0.25", "2", "16", "64"))
#' essential_agreement(x, y)
#' # TRUE FALSE FALSE TRUE
essential_agreement <- function(x,
                                y,
                                coerce_mic = TRUE,
                                mode = "categorical") {
  if (any(!AMR::is.mic(c(x, y))) & !coerce_mic) {
    stop("Both MIC inputs to essential_agreement must be AMR::mic.
Convert using AMR::as.mic() with or without molMIC::force_mic().")
  }

  # if (any(!AMR::is.mic(c(x, y)))) {
  #   x <- AMR::as.mic(x)
  #   y <- AMR::as.mic(y)
  # }

  if (mode == "categorical") {

    x <- force_mic(mic_uncensor(x))
    y <- force_mic(mic_uncensor(y))

    index_diff <- purrr::map2_lgl(x,
                                  y,
                                  \(.x, .y) {
                                    range <- mic_range()
                                    range <- c(max(range) * 2,
                                               range,
                                               min(range) / 2)
                                    if (any(is.na(c(.x, .y)))) {
                                      return(NA)
                                    }
                                    if(abs(which(range == .x) - which(range == .y)) > 1) {
                                      return(FALSE)
                                    }
                                    return(TRUE)
                                  })

    return(index_diff)

  }
  if (mode == "numerical") {
    x <- as.numeric(mic_uncensor(x))
    y <- as.numeric(mic_uncensor(y))

    frac <- x / y
    return(as.logical(
      dplyr::case_when(
        is.na(frac) ~ NA,
        frac > 2.0 ~ FALSE,
        frac < 0.5 ~ FALSE,
        frac != 1.0 ~ TRUE,
        TRUE ~ TRUE)))
  }
  stop("Mode must be categorical or numerical")
}

#' Perform an MIC validation experiment
#'
#' @param gold_standard vector of MICs to compare against.
#' @param test vector of MICs that are under investigation
#' @param ab character vector (same length as MIC) of antibiotic names (optional)
#' @param mo character vector (same length as MIC) of microorganism names (optional)
#' @param simplify if TRUE, MIC values will be coerced into the closest halving
#' dilution (e.g., 0.55 will be converted to 0.5)
#'
#' @return S3 mic_validation object
#' @export
compare_mic <- function(gold_standard,
                        test,
                        ab = NULL,
                        mo = NULL,
                        simplify = TRUE) {
  gold_standard_mod <- gold_standard |>
    force_mic(levels_from_AMR = !simplify) |>
    AMR::as.mic()

  test_mod <- test |>
    force_mic(levels_from_AMR = !simplify) |>
    AMR::as.mic()

  output <- data.frame(
    gold_standard = gold_standard_mod,
    test = test_mod,
    essential_agreement = essential_agreement(gold_standard_mod, test_mod)
  )

  if (!is.null(ab)) {
    output["ab"] <- ab
  }

  if (!is.null(mo)) {
    output["mo"] <- mo
  }

  if (!is.null(ab) & !is.null(mo)) {
    gold_standard_sir <- purrr::pmap_vec(list(gold_standard_mod,
                                     mo,
                                     ab), \(x, mo, ab) AMR::as.sir(x, mo = mo, ab = ab))
    test_sir <- purrr::pmap_vec(list(test_mod,
                                     mo,
                                     ab), \(x, mo, ab) AMR::as.sir(x, mo = mo, ab = ab))
    output[["gold_standard_sir"]] <- gold_standard_sir
    output[["test_sir"]] <- test_sir
    output["error"] <- compare_sir(gold_standard_sir,
                                   test_sir)
  }

  class(output) <- append(class(output), "mic_validation", 0)
  output
}

plot_mic_validation_single_ab <- function(x, match_axes, ...) {
  x <- x |>
    dplyr::group_by(.data[["gold_standard"]],
                    .data[["test"]],
                    .data[["essential_agreement"]]) |>
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
    ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(fill=NA))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::xlab("Gold standard MIC (mg/L)") +
    ggplot2::ylab("Test (mg/L)")

  if (match_axes) {
    x <- x + ggplot2::scale_x_discrete(drop = FALSE)
    x <- x + ggplot2::scale_y_discrete(drop = FALSE)
  }
  x
}

plot_mic_validation_multi_ab <- function(x, match_axes, ...) {
  x <- x |>
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
    #ggplot2::theme(legend.position = c(0.9,0.2)) +
    ggplot2::xlab("Gold standard MIC (mg/L)") +
    ggplot2::ylab("Test (mg/L)")

  if (match_axes) {
    x <- x + ggplot2::scale_x_discrete(drop = FALSE)
    x <- x + ggplot2::scale_y_discrete(drop = FALSE)
  }
  x
}

#' Plot MIC validation results
#'
#' @param x object generated using compare_mic
#' @param match_axes Same x and y axis
#' @param ... additional arguments
#'
#' @export
plot.mic_validation <- function(x, match_axes = TRUE, ...) {

  if (match_axes) {
    all_levels <- levels(AMR::as.mic(NA))

    keep_these <- levels(droplevels(x[["test"]]))
    drop_these <- all_levels[!all_levels %in% keep_these]

    x[["gold_standard"]] <- forcats::fct_drop(x[["gold_standard"]], only = drop_these)


    keep_these <- levels(droplevels(x[["gold_standard"]]))
    drop_these <- all_levels[!all_levels %in% keep_these]

    x[["test"]] <- forcats::fct_drop(x[["test"]], only = drop_these)

    class(x[["test"]]) <- append(class(x[["test"]]), values = "mic", after = 0)
    class(x[["gold_standard"]]) <- append(class(x[["gold_standard"]]), values = "mic", after = 0)
  }

  if (!"ab" %in% colnames(x) | length(unique(x[["ab"]])) == 1) {
    plot_mic_validation_single_ab(x, match_axes, ...)
  }
  else {
    plot_mic_validation_multi_ab(x, match_axes, ...)
  }
}

#' Summary of MIC validation results
#'
#' @param object S3 mic_validation object
#' @param ... further optional parameters
#'
#' @export
summary.mic_validation <- function(object, ...) {
  if (!"ab" %in% colnames(object) & !"mo" %in% colnames(object)) {
    return(
      list(EA = sum(object$essential_agreement == TRUE) / length(object$essential_agreement),
           bias = bias(object$gold_standard, object$test))
    )
  }

  if ("ab" %in% colnames(object) & !"mo" %in% colnames(object)) {
    return(
      object |>
        dplyr::group_by(.data[["ab"]]) |>
        dplyr::summarise(EA = sum(.data[["essential_agreement"]] == TRUE) / length(.data[["essential_agreement"]]),
                         bias = bias(object$gold_standard, object$test))
    )
  }

  if ("ab" %in% colnames(object) & "mo" %in% colnames(object)) {
    return(
      object |>
        dplyr::group_by(.data[["ab"]], .data[["mo"]]) |>
        dplyr::summarise(EA = sum(.data[["essential_agreement"]] == TRUE) / length(.data[["essential_agreement"]]),
                         `minor error (%)` = sum(.data[["error"]] == "m") / length(.data[["error"]]) * 100,
                         `major error (%)` = sum(.data[["error"]] == "M") / length(.data[["error"]]) * 100,
                         `very major error (%)` = sum(.data[["error"]] == "vM") /length(.data[["error"]]) * 100,
                         `minor error (n)` = sum(.data[["error"]] == "m"),
                         `major error (n)` = sum(.data[["error"]] == "M"),
                         `very major error (n)` = sum(.data[["error"]] == "vM"))
    )
  }
}

#' Check that MIC is within QC range
#'
#' @param measurement measured QC MIC
#' @param strain control strain identifier (usually ATCC)
#' @param ab antibiotic name (will be coerced to AMR::as.ab)
#' @param ignore_na ignores NA (returns TRUE)
#' @param guideline Guideline to use (EUCAST or CLSI)
#' @param year Guideline year (version)
#' @return logical vector
#' @export
#'
#' @examples
#' qc_in_range(AMR::as.mic(0.5), 25922, "GEN") == TRUE
#' qc_in_range(AMR::as.mic(8.0), 25922, "GEN") == FALSE
qc_in_range <- Vectorize(
    function(measurement,
             strain,
             ab,
             ignore_na = TRUE,
             guideline = "EUCAST",
             year = "2023") {
    if (is.na(measurement) | is.na(strain) | is.na(ab)) {
      if (ignore_na) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    if (!inherits(measurement, "mic") & !inherits(measurement, "disk")) {
      stop("Measurement must be AMR::mic or AMR::disk")
    }
    ab <- as.character(AMR::as.ab(ab))

    strain <- tolower(strain)
    if (!startsWith(strain, "atcc")) {
      strain <- paste0("atcc", strain)
    }
    qc_subset <- QC_table[QC_table$STRAIN == strain & QC_table$ANTIBIOTIC == ab & QC_table$GUIDELINE == guideline & QC_table$YEAR == year & QC_table$METHOD == "MIC", ]

    if (nrow(qc_subset) == 0) {
      # no QC info in table
      return(NA)
    }

    stopifnot(nrow(qc_subset) == 1)

    if (inherits(measurement, "mic")) {
      measurement <- mic_uncensor(measurement)

      if (measurement < AMR::as.mic(qc_subset[, "MINIMUM"])) {
        return(FALSE)
      }
      if (measurement > AMR::as.mic(qc_subset[, "MAXIMUM"])) {
        return(FALSE)
      }
      return(TRUE)
    }
    },
    vectorize.args = c("measurement", "strain", "ab")
)

#' Check that QC measurement is at the required target
#'
#' @param measurement measured QC MIC
#' @param strain control strain identifier (usually ATCC)
#' @param ab antibiotic name (will be coerced to AMR::as.ab)
#' @param ignore_na ignores NA (returns TRUE)
#' @param guideline Guideline to use (EUCAST or CLSI)
#' @param year Guideline year (version)
#'
#' @return logical vector
#' @export
#'
#' @examples
#' qc_on_target(AMR::as.mic(0.5), 25922, "GEN") == TRUE
qc_on_target <- Vectorize(
  function(measurement,
           strain,
           ab,
           ignore_na = TRUE,
           guideline = "EUCAST",
           year = "2023") {
    if (is.na(measurement) | is.na(strain) | is.na(ab)) {
      if (ignore_na) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    if (!inherits(measurement, "mic") & !inherits(measurement, "disk")) {
      stop("Measurement must be AMR::mic or AMR::disk")
    }
    ab <- as.character(AMR::as.ab(ab))

    strain <- tolower(strain)
    if (!startsWith(strain, "atcc")) {
      strain <- paste0("atcc", strain)
    }

    qc_subset <- QC_table[QC_table$STRAIN == strain & QC_table$ANTIBIOTIC == ab & QC_table$GUIDELINE == guideline & QC_table$YEAR == year & QC_table$METHOD == "MIC", ]

    if (nrow(qc_subset) == 0) {
      # no QC info in table
      return(NA)
    }

    stopifnot(nrow(qc_subset) == 1)

    if (inherits(measurement, "mic")) {
      measurement <- mic_uncensor(measurement)

      if (measurement < qc_subset[, "MINIMUM_TARGET"]) {
        return(FALSE)
      }
      if (measurement > qc_subset[, "MAXIMUM_TARGET"]) {
        return(FALSE)
      }
      return(TRUE)
    }
  },
  vectorize.args = c("measurement", "strain", "ab")
)

#' Standardise MIC to control strain
#'
#' @param test_measurement Measured MIC to standardise
#' @param qc_measurement Measured QC MIC to standardise to
#' @param strain control strain identifier (usually ATCC)
#' @param ab antibiotic name (will be coerced to AMR::as.ab)
#' @param prefer_upper Where the target MIC is a range, prefer the upper value
#' in the range
#' @param ignore_na Ignore NA (returns AMR::NA_mic_)
#' @param guideline Guideline to use (EUCAST or CLSI)
#' @param year Guideline year (version)
#' @param force Force into MIC-compatible format after standardisation
#'
#' @return AMR::mic vector
#' @export
#'
#' @examples
#' standardise_mic(
#'   test_measurement = c(AMR::as.mic(">8.0"),
#'                        AMR::as.mic(4.0),
#'                        AMR::as.mic(2)),
#'   qc_measurement = c(AMR::as.mic(1),
#'                      AMR::as.mic(0.5),
#'                      AMR::as.mic(1)),
#'   strain = 25922,
#'   ab = AMR::as.ab("GEN"))
standardise_mic <- function(test_measurement,
                            qc_measurement,
                            strain,
                            ab,
                            prefer_upper = FALSE,
                            ignore_na = TRUE,
                            guideline = "EUCAST",
                            year = "2023",
                            force = TRUE) {
  standardise_mic_vectorised <- Vectorize(
    function(test_measurement,
             qc_measurement,
             strain,
             ab,
             prefer_upper,
             ignore_na,
             guideline,
             year,
             force) {
      if (is.na(test_measurement) | is.na(qc_measurement) | is.na(strain) | is.na(ab)) {
        if (ignore_na) {
          return(NA)
        } else {
          return(NA)
        }
      }
      if (!inherits(test_measurement, "mic") | !inherits(qc_measurement, "mic")) {
        stop("Measurements must be AMR::mic")
      }
      ab <- as.character(AMR::as.ab(ab))

      mic_char <- as.character(test_measurement)
      if (grepl("<|>", mic_char)) {
        return(as.character(test_measurement))
      }

      if (qc_on_target(qc_measurement, strain, ab)) {
        return(as.character(test_measurement))
      }

      if (!qc_in_range(qc_measurement, strain, ab)) {
        warning("QC not in range, standardise_mic may not be appropriate")
      }

      strain <- tolower(strain)
      if (!startsWith(strain, "atcc")) {
        strain <- paste0("atcc", strain)
      }
      qc_subset <- QC_table[QC_table$STRAIN == strain & QC_table$ANTIBIOTIC == ab & QC_table$GUIDELINE == guideline & QC_table$YEAR == year & QC_table$METHOD == "MIC", ]

      if (nrow(qc_subset) == 0) {
        # no QC info in table
        return(NA)
      }

      stopifnot(nrow(qc_subset) == 1)

      if (qc_subset[1, "TARGET_ESTIMATED"]) {
        warning(paste("Target for", ab, "is estimated as mid-point between upper
and lower range. EUCAST/CLSI guidance may be different."))
      }

      scale_factor_lower <- qc_subset[, "MINIMUM_TARGET"]
      scale_factor_upper <- qc_subset[, "MAXIMUM_TARGET"]
      if (scale_factor_lower == scale_factor_upper) {
        scale_factor <- scale_factor_lower
      } else{
        if (prefer_upper) {
          scale_factor <- scale_factor_upper
        } else {
          scale_factor <- scale_factor_lower
        }
      }

      as.character(test_measurement * (scale_factor / qc_measurement))
    },
    vectorize.args = c("test_measurement",
                       "qc_measurement",
                       "strain",
                       "ab")
  )
  if (force) {
    return(
      AMR::as.mic(force_mic(standardise_mic_vectorised(test_measurement,
                                                       qc_measurement,
                                                       strain,
                                                       ab,
                                                       prefer_upper,
                                                       ignore_na,
                                                       guideline,
                                                       year)
      )
    ))
  }
  AMR::as.mic(standardise_mic_vectorised(test_measurement,
                                         qc_measurement,
                                         strain,
                                         ab,
                                         prefer_upper,
                                         ignore_na,
                                         guideline,
                                         year))
}

compare_sir <- function(gold_standard, test) {
  if (length(gold_standard) != length(test)) {
    stop("Inputs to compare_sir must be same length")
  }
  return(
    factor(
      dplyr::case_when(
        gold_standard == "S" & test == "I" ~ "m",
        gold_standard == "R" & test == "I" ~ "m",
        gold_standard == "I" & test == "S" ~ "m",
        gold_standard == "I" & test == "R" ~ "m",
        gold_standard == "S" & test == "R" ~ "M",
        gold_standard == "R" & test == "S" ~ "vM",
        TRUE ~ NA
      )
    )
  )
}

bias <- function(gold_standard, test) {
  if (length(gold_standard) != length(test)) {
    stop("Inputs to bias must be same length")
  }
  bias_above <- sum(test > gold_standard)
  bias_below <- sum(test < gold_standard)

  n <- length(gold_standard)
  return(bias_above / n * 100 - bias_below / n * 100)
}
