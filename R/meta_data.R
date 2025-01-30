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
#' @description
#' This function helps extract MICs from a database of results. It is compatible
#' with the PATRIC meta data format when used on a tidy_patric_db object,
#' created using tidy_patric_db().
#'
#' If more than one MIC is present for a particular observation, the function
#' can return the higher MIC by setting prefer_high_mic = TRUE. If
#' prefer_high_mic = FALSE, the lower MIC will be returned.
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
#'
#' @param mic vector of MICs to uncensor; will be coerced to MIC using AMR::as.mic
#' @param method method to uncensor MICs (scale, simple, or bootstrap)
#' @param scale scalar to multiply or divide MIC by (for method = scale)
#' @param ab antibiotic name (for method = bootstrap)
#' @param mo microorganism name (for method = bootstrap)
#' @param distros dataframe of epidemiological distributions (only used,
#' optionally, for method = bootstrap)
#'
#' @return vector of MICs in AMR::mic format
#'
#' @details
#' Censored MIC data is generally unsuitable for modelling without some
#' conversion of censored data. The default behaviour (method = scale) is to
#' halve MICs under the limit of detection (<=) and double MICs above the limit
#' of detection (>). When used with method = simple, this function effectively
#' just removes the censoring symbols, e.g., <=2 becomes 2, and >64 becomes 64.
#'
#' The bootstrap method is the more complex of the three available methods. It
#' attempts to use a second (uncensored) MIC distribution to sample values in
#' the censored range. These values are then used to populate and uncensor
#' the MIC data provided as input (mic). The second (uncensored) MIC
#' distribution is ideally provided from similar experimental conditions.
#' Alternatively, epidemiological distributions can be used. These distributions
#' should be provided as a dataframe to the distros argument. The format for
#' this dataframe is inspired by the EUCAST epidemiological distributions, see:
#' https://www.eucast.org/mic_and_zone_distributions_and_ecoffs. The dataframe
#' should contain columns for antimicrobial (converted using AMR::as.ab),
#' organism (converted using AMR::as.mo), and MIC concentrations. An example
#' is provided in the 'ecoffs' dataset available with this pacakge. Currently,
#' only Escherichia coli is available in this dataset. Each observation (row)
#' consists of the frequency a particular MIC concentration is observed in the
#' distribution. If such a dataframe is not provided to distros, the function
#' will attempt to use 'ecoffs', but remains limited to E. coli.
#'
#' @references https://www.eucast.org/mic_and_zone_distributions_and_ecoffs
#'
#' @export
#'
#' @examples
#' mic_uncensor(c(">64.0", "<0.25", "8.0"), method = "scale", scale = 2)
mic_uncensor <- function(mic,
                         method = "scale",
                         scale = 2,
                         ab = NULL,
                         mo = NULL,
                         distros = NULL) {
  if (method == "scale") {
    return(mic_uncensor_scale(mic, scale))
  }
  if (method == "simple") {
    return(mic_uncensor_simple(mic))
  }
  if (method == "bootstrap") {
    if (any(is.null(ab),
            is.null(mo))) {
      stop("For method = bootstrap, both ab and mo must be provided")
    }
    return(AMR::as.mic(mic_uncensor_bootstrap(mic = mic, ab = ab, mo = mo,
                                              distros = distros)))
  }
  stop("Method must be scale, simple or bootstrap")
}

mic_uncensor_simple <- function(mic) {
  AMR::as.mic(as.numeric(AMR::as.mic(mic)))
}

mic_uncensor_scale <- function(mic, scale) {
  suppressWarnings(
    dplyr::case_when(
      stringr::str_detect(mic, ">") ~ AMR::as.mic(AMR::as.mic(stringr::str_remove(force_mic(mic), ">")) * scale),
      stringr::str_detect(mic, "<") ~ AMR::as.mic(AMR::as.mic(stringr::str_remove(force_mic(mic), "<")) / scale),
      .default = AMR::as.mic(force_mic(mic)))
  )
}

mic_uncensor_bootstrap <- Vectorize(function(mic, ab, mo, distros) {
  if (is.na(mic)) {
    return(AMR::NA_mic_)
  }
  mic <- gsub("=", "", mic)
  if (!startsWith(mic, "<") & !startsWith(mic, ">")) {
    return(AMR::as.mic(mic))
  }

  ab <- AMR::as.ab(ab)
  mic <- AMR::as.mic(mic)
  if (!is.null(mo)) {
    mo <- AMR::as.mo(mo)
  }

  if (is.null(distros)) {
    message(
    "No custom distributions provided to mic_uncensor using bootstrap method.
Using built-in epidemiological distributions from EUCAST:
https://www.eucast.org/mic_and_zone_distributions_and_ecoffs.
Note: Only Escherichia coli is currently supported.")
    distros <- ecoffs
  }

  ecoff_sub <- distros[distros['antibiotic'] == ab & distros['organism'] == mo,]
  ecoff_mics <- ecoff_sub[as.character(mic_range())]

  if (nrow(ecoff_sub) == 0) {
    warning("No ECOFFs available for ", ab)
    return(mic)
  }

  if (startsWith(as.character(mic), "<")) {
    ecoff_mics[AMR::as.mic(names(ecoff_mics)) > mic] <- 0
  }
  if (startsWith(as.character(mic), ">")) {
    ecoff_mics[AMR::as.mic(names(ecoff_mics)) <= mic] <- 0
  }

  ecoff_mics <- sapply(ecoff_mics, as.numeric)

  if (sum(ecoff_mics) == 0) {
    warning(ab, "appears to be below or above epidemiological distribution,
performing a simple mic uncensor.")
    return(mic_uncensor_simple(mic))
  }

  return(AMR::as.mic(sample(names(ecoff_mics), size = 1, replace=T, prob = ecoff_mics)))
},
USE.NAMES = F,
vectorize.args = c("mic", "ab", "mo"))

#' Censor MIC values
#'
#' @param mic MIC (coercible to AMR::as.mic)
#' @param ab antibiotic name (coercible to AMR::as.ab)
#' @param mo microorganism name (coercible to AMR::as.mo)
#' @param rules censor rules - named list of pathogen (in AMR::as.mo code) to
#' antibiotic (in AMR::as.ab code) to censoring rules. The censoring rules
#' should provide a min or max value to censor MICs to. See example for more.
#'
#' @return censored MIC values (S3 mic class)
#'
#' @description
#' MIC datasets often arise from different laboratories or experimental
#' conditions. In practice, this means that there can be different levels of
#' censoring (<= and >) within the data. This function can be used to harmonise
#' the dataset to a single level of censoring. The function requires a set of
#' rules that specify the censoring levels (see example).
#'
#' @export
#'
#' @examples
#' example_rules <- list("B_ESCHR_COLI" = list(
#'   "AMK" = list(min = 2, max = 32),
#'   "CHL" = list(min = 4, max = 64),
#'   "GEN" = list(min = 1, max = 16),
#'   "CIP" = list(min = 0.015, max = 4),
#'   "MEM" = list(min = 0.016, max = 16),
#'   "AMX" = list(min = 2, max = 64),
#'   "AMC" = list(min = 2, max = 64),
#'   "FEP" = list(min = 0.5, max = 64),
#'   "CAZ" = list(min = 1, max = 128),
#'   "TGC" = list(min = 0.25, max = 1)
#'   ))
#'
#' mic_censor(AMR::as.mic(512),
#'            "AMK",
#'            "B_ESCHR_COLI",
#'            example_rules) == AMR::as.mic(">32")
mic_censor <- function(mic, ab, mo, rules) {
  mic_censor_vectorize <- Vectorize(
    function(mic, ab, mo, rules) {
      mic <- AMR::as.mic(mic)
      min_thresh <- rules[[mo]][[ab]][["min"]]
      min_thresh <- ifelse(is.null(min_thresh),
                           -Inf,
                           min_thresh)
      max_thresh <- rules[[mo]][[ab]][["max"]]
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
    USE.NAMES = FALSE,
    vectorize.args = c("mic", "ab", "mo")
  )

  AMR::as.mic(mic_censor_vectorize(mic, ab, mo, rules))
}

#' Generate dilution series
#'
#' @param start starting (highest) concentration
#' @param dilutions number of dilutions
#' @param min minimum (lowest) concentration
#' @param precise force range to be high precision (not usually desired
#' behaviour)
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
#' @description
#' Convert a value that is "almost" an MIC into a valid MIC value.
#'
#' @param value vector of MIC-like values (numeric or character)
#' @param levels_from_AMR conform to AMR::as.mic levels
#' @param max_conc maximum concentration to force to
#' @param min_conc minimum concentration to force to
#' @param method method to use when forcing MICs (closest or round_up)
#' @param prefer where value is in between MIC (e.g., 24mg/L) chose the higher
#' MIC ("max") or lower MIC ("min"); only applies to method = "closest"
#'
#' @details
#' Some experimental or analytical conditions measure MIC (or surrogate) in a
#' way that does not fully conform to traditional MIC levels
#' (i.e., concentrations). This function allows these values to be coerced into
#' an MIC value that is compatible with the AMR::mic class. When using method =
#' "closest", the function will choose the closest MIC value to the input value
#' (e.g., 2.45 will be coerced to 2). When using method = "round up", the
#' function will round up to the next highest MIC value (e.g., 2.45 will be
#' coerced to 4). "Round up" is technically the correct approach if the input
#' value was generated from an experiment that censored between concentrations
#' (e.g., broth or agar dilution). However, "closest" may be more appropriate in
#' some cases.
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
                      method = "closest",
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

    if (!method %in% c("closest", "round up")) {
      stop("Method must be closest or round up")
    }

    if (method == "closest") {
      positions <- which(abs(mic_vector - inner_val) == min(abs(mic_vector - inner_val)))
    }

    if (method == "round up") {
      if (inner_val %in% mic_vector) {
        positions <- which(mic_vector == inner_val)
      } else {
        positions <- min(which(inner_val < mic_vector))
      }
    }

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
#' @description
#' Essential agreement calculation for comparing two MIC vectors.
#'
#' @param x AMR::mic or coercible
#' @param y AMR::mic or coercible
#' @param coerce_mic convert to AMR::mic
#' @param mode Categorical or numeric
#' @return logical vector
#' @export
#'
#' @details
#' Essential agreement is a central concept in the comparison of two sets of MIC
#' values. It is most often used when validating a new method against a gold
#' standard. This function reliably performs essential agreement in line with
#' ISO 20776-2:2021. The function can be used in two modes: categorical and
#' numeric. In categorical mode, the function will use traditional MIC
#' concentrations to determine the MIC (therefore it will use force_mic() to
#' convert both x and y to a clean MIC -- see ?force_mic()). In numeric mode,
#' the function will compare the ratio of the two MICs. In most cases,
#' categorical mode provides more reliable results. Values within +/- 2
#' dilutions are considered to be in essential agreement.
#'
#' @references
#' International Organization for Standardization. ISO 20776-2:2021
#' Available from: https://www.iso.org/standard/79377.html
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
Convert using AMR::as.mic() with or without MIC::force_mic().")
  }

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

#' Compare and validate MIC values
#'
#' @param gold_standard vector of MICs to compare against.
#' @param test vector of MICs that are under investigation
#' @param ab character vector (same length as MIC) of antibiotic names (optional)
#' @param mo character vector (same length as MIC) of microorganism names (optional)
#' @param accept_ecoff if TRUE, ECOFFs will be used when no clinical breakpoints are available
#' @param simplify if TRUE, MIC values will be coerced into the closest halving
#' dilution (e.g., 0.55 will be converted to 0.5)
#'
#' @return S3 mic_validation object
#'
#' @description
#' This function compares an vector of MIC values to another. Generally, this is
#' in the context of a validation experiment -- an investigational assay or
#' method (the "test") is compared to a gold standard. The rules used by this
#' function are in line with "ISO 20776-2:2021 Part 2: Evaluation of performance
#' of antimicrobial susceptibility test devices against reference broth
#' micro-dilution."
#'
#' There are two levels of detail that are provided. If only the MIC values are
#' provided, the function will look for essential agreement between the two sets
#' of MIC. If the organism and antibiotic arguments are provided, the function
#' will also calculate the categorical agreement using EUCAST breakpoints (or,
#' if breakpoint not available and accept_ecoff = TRUE, ECOFFs).
#'
#' The function returns a special dataframe of results, which is also an
#' mic_validation object. This object can be summarised using summary() for
#' summary metrics, plotted using plot() for an essential agreement confusion
#' matrix, and tabulated using table().
#'
#' @export
#'
#' @examples
#' # Just using MIC values only
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' summary(val)
#'
#' # Using MIC values and antibiotic and organism names
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' ab <- c("AMK", "AMK", "AMK", "AMK")
#' mo <- c("B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI")
#' val <- compare_mic(gold_standard, test, ab, mo)
#' "error" %in% names(val)  # val now has categorical agreement
compare_mic <- function(gold_standard,
                        test,
                        ab = NULL,
                        mo = NULL,
                        accept_ecoff = FALSE,
                        simplify = TRUE) {
  if (length(gold_standard) != length(test)) {
    stop("Gold standard and test must be the same length")
  }

  gold_standard_mod <- gold_standard |>
    force_mic(levels_from_AMR = !simplify) |>
    AMR::as.mic()

  test_mod <- test |>
    force_mic(levels_from_AMR = !simplify) |>
    AMR::as.mic()

  output <- list(
    gold_standard = gold_standard_mod,
    test = test_mod,
    essential_agreement = factor(essential_agreement(gold_standard_mod, test_mod),
                                 levels = c(FALSE, TRUE))
  )

  if (!is.null(ab)) {
    if (length(ab) == 1) {
      ab <- rep(ab, length(gold_standard))
    }
    if (length(ab) > 1 & length(ab) != length(gold_standard)) {
      stop("Antibiotic names must be the same length as MIC values, or single value")
    }
    output[["ab"]] <- ab
  }

  if (!is.null(mo)) {
    if (length(mo) == 1) {
      mo <- rep(mo, length(gold_standard))
    }
    if (length(mo) > 1 & length(mo) != length(gold_standard)) {
      stop("Microorganism names must be the same length as MIC values, or single value")
    }
    output[["mo"]] <- mo
  }

  if (!is.null(ab) & !is.null(mo)) {
    gold_standard_sir <- purrr::pmap_vec(list(gold_standard_mod,
                                              mo,
                                              ab), \(x, mo, ab) AMR::as.sir(x, mo = mo, ab = ab))
    gold_standard_sir <- purrr::pmap_vec(list(gold_standard_mod,
                                              mo,
                                              ab,
                                              gold_standard_sir), \(x, mo, ab, sir) {
                                                if (is.na(sir)) {
                                                  if (accept_ecoff) {
                                                    message(paste("No clinical breakpoint for:", ab, mo, "using ECOFF"))
                                                    AMR::as.sir(x, mo = mo, ab = ab,
                                                                breakpoint_type = "ECOFF")
                                                  } else {
                                                    message(paste("No clinical breakpoint for:", ab, mo, "(use accept_ecoff to use ECOFF)"))
                                                    sir
                                                  }
                                                } else {
                                                  sir
                                                }
                                              })
    test_sir <- purrr::pmap_vec(list(test_mod,
                                     mo,
                                     ab), \(x, mo, ab) AMR::as.sir(x, mo = mo, ab = ab))
    test_sir <- purrr::pmap_vec(list(test_mod,
                                     mo,
                                     ab,
                                     test_sir), \(x, mo, ab, sir) {
                                       if (is.na(sir)) {
                                         if (accept_ecoff) {
                                           message(paste("No clinical breakpoint for:", ab, mo, "using ECOFF"))
                                           AMR::as.sir(x, mo = mo, ab = ab,
                                                       breakpoint_type = "ECOFF")
                                         } else {
                                           message(paste("No clinical breakpoint for:", ab, mo, "(use accept_ecoff to use ECOFF)"))
                                           sir
                                         }
                                       } else {
                                         sir
                                       }
                                     })
    output[["gold_standard_sir"]] <- gold_standard_sir
    output[["test_sir"]] <- test_sir
    output[["error"]] <- compare_sir(gold_standard_sir,
                                   test_sir)
  }

  class(output) <- append(class(output), "mic_validation", 0)
  output
}

#' Print MIC validation object
#'
#' @param x mic_validation object
#' @param ... additional arguments
#'
#' @return character
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' print(val)
print.mic_validation <- function(x, ...) {
  if (!is.null(x$ab) & !is.null(x$mo)) {
    antibiotic_str <- paste(unique(x$ab), collapse = ", ")
    organism_str <- paste(unique(x$mo), collapse = ", ")
    return(
      cat(glue::glue("MIC validation object with {length(x$gold_standard)} observations
                     Agreement type: essential and categorical
                     Antibiotics: {antibiotic_str}
                     Organisms: {organism_str}"))
      )
    }

  return(
    cat(glue::glue("MIC validation object with {length(x$gold_standard)} observations
                   Agreement type: essential"))
  )

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
    ggplot2::geom_tile(alpha=1, show.legend = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label=.data[["n"]]), show.legend = TRUE) +
    ggplot2::scale_fill_gradient(low="white", high="#009194") +
    ggplot2::scale_fill_manual(values=c("red", "black"), aesthetics = "color", drop = FALSE) +
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

plot_mic_validation_multi_ab <- function(x, match_axes,
                                         facet_wrap_ncol,
                                         facet_wrap_nrow,
                                         ...) {
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
    ggplot2::scale_fill_manual(values=c("red", "black"), aesthetics = "color", drop = FALSE)

    if (any(!is.null(c(facet_wrap_ncol, facet_wrap_nrow)))) {
      x <- x + lemon::facet_rep_wrap(~ .data[["ab"]],
                                     nrow = facet_wrap_nrow,
                                     ncol = facet_wrap_ncol,
                                     repeat.tick.labels = TRUE)
    }

    x <- x +
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

#' S breakpoint for MIC
#'
#' @param mo mo name (coerced using AMR::as.mo)
#' @param ab ab name (coerced using AMR::as.ab)
#' @param accept_ecoff if TRUE, ECOFFs will be used when no clinical breakpoints are available
#' @param ... additional arguments to pass to AMR::as.sir, which is used to
#' calculate the S breakpoint
#'
#' @return MIC value
#' @export
#'
#' @examples
#' mic_s_breakpoint("B_ESCHR_COLI", "AMK")
#' mic_s_breakpoint("B_ESCHR_COLI", "CHL", accept_ecoff = TRUE)
mic_s_breakpoint <- function(mo, ab, accept_ecoff = FALSE, ...) {
  mic_range <- mic_range()
  sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                           mo = mo,
                           ab = ab,
                           ...)
  if (any(is.na(sir_range))) {
    if (accept_ecoff) {
      sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                               mo = mo,
                               ab = ab,
                               breakpoint_type = "ECOFF",
                               ...)
    } else {
      stop("No clinical breakpoints available for ", ab, mo, "(consider using accept_ecoff if appropriate.")
    }
  }
  for (i in seq_along(sir_range)) {
    if (sir_range[i] == "S") {
      return(AMR::as.mic(mic_range[i]))
    }
  }
}

#' R breakpoint for MIC
#'
#' @param mo mo name (coerced using AMR::as.mo)
#' @param ab ab name (coerced using AMR::as.ab)
#' @param accept_ecoff if TRUE, ECOFFs will be used when no clinical breakpoints are available
#' @param ... additional arguments to pass to AMR::as.sir, which is used to
#' calculate the R breakpoint
#'
#' @return MIC value
#' @export
#'
#' @examples
#' mic_r_breakpoint("B_ESCHR_COLI", "AMK")
#' mic_r_breakpoint("B_ESCHR_COLI", "CHL", accept_ecoff = TRUE)
mic_r_breakpoint <- function(mo, ab, accept_ecoff = FALSE, ...) {
  mic_range <- rev(mic_range())
  sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                           mo = mo,
                           ab = ab,
                           ...)
  if (any(is.na(sir_range))) {
    if (accept_ecoff) {
      sir_range <- AMR::as.sir(AMR::as.mic(mic_range),
                               mo = mo,
                               ab = ab,
                               breakpoint_type = "ECOFF",
                               ...)
    } else {
      stop("No clinical breakpoints available for ", ab, mo, "(consider using accept_ecoff if appropriate.")
    }
  }
  for (i in seq_along(sir_range)) {
    if (sir_range[i] == "R") {
      return(AMR::as.mic(mic_range[i-1]))
    }
  }
}

match_levels <- function(x, match_to) {
  if (!AMR::is.mic(x) | !AMR::is.mic(match_to)) {
    stop("Both x and match_to must be AMR::mic objects")
  }
  all_levels <- levels(AMR::as.mic(NA))

  keep_these <- levels(droplevels(match_to))
  drop_these <- all_levels[!all_levels %in% keep_these]

  x <- forcats::fct_drop(x, only = drop_these)
  class(x) <- append(class(x), "mic", after = 0)
  x
}

#' Fill MIC dilution levels
#'
#' @param x MIC vector
#' @param cap_upper If True, will the top level will be the highest MIC dilution in x
#' @param cap_lower If True, will the bottom level will be the lowest MIC dilution in x
#' @param as.mic By default, returns an ordered factor. Set as.mic = TRUE to return as AMR::mic
#'
#' @return ordered factor (or AMR::mic if as.mic = TRUE)
#' @export
#'
#' @examples
#' # use in combination with droplevels to clean up levels:
#' x <- AMR::as.mic(c("<0.25", "8", "64", ">64"))
#' x <- droplevels(x)
#' fill_dilution_levels(x)
fill_dilution_levels <- function(x,
                                 cap_upper = TRUE,
                                 cap_lower = TRUE,
                                 as.mic = TRUE) {
  mic_conc_range <- AMR::as.mic(mic_range())
  if (cap_lower) {
    mic_conc_range <- mic_conc_range[mic_conc_range >= min(AMR::as.mic(x))]
  }
  if (cap_upper) {
    mic_conc_range <- mic_conc_range[mic_conc_range <= max(AMR::as.mic(x))]
  }
  x <- forcats::fct_expand(x, as.character(mic_conc_range))
  new_levels <- levels(x)[order(match(levels(x),
                                      as.character(levels(AMR::as.mic(NA)))))]
  x <- ordered(x, levels = new_levels)
  if (as.mic) {
    class(x) <- append(class(x), "mic", after = 0)
  }
  x
}

#' Plot MIC validation results
#'
#' @param x object generated using compare_mic
#' @param match_axes Same x and y axis
#' @param add_missing_dilutions Axes will include dilutions that are not
#' @param facet_wrap_ncol Facet wrap into n columns by antimicrobial (optional,
#' only available when more than one antimicrobial in validation)
#' @param facet_wrap_nrow Facet wrap into n rows by antimicrobial (optional,
#' only available when more than one antimicrobial in validation)
#' represented in the data, based on a series of dilutions generated using mic_range().
#' @param ... additional arguments
#'
#' @return ggplot object
#'
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' plot(val)
#'
#' # works with validation that includes categorical agreement
#' # categorical agreement is ignored
#' ab <- c("AMK", "AMK", "AMK", "AMK")
#' mo <- c("B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI", "B_ESCHR_COLI")
#' val <- compare_mic(gold_standard, test, ab, mo)
#' plot(val)
#'
#' # if the validation contains multiple antibiotics, i.e.,
#' ab <- c("CIP", "CIP", "AMK", "AMK")
#' val <- compare_mic(gold_standard, test, ab, mo)
#' # the following will plot all antibiotics in a single plot (pooled results)
#' plot(val)
#' # use the faceting arguments to split the plot by antibiotic
#' plot(val, facet_wrap_ncol = 2)
plot.mic_validation <- function(x,
                                match_axes = TRUE,
                                add_missing_dilutions = TRUE,
                                facet_wrap_ncol = NULL,
                                facet_wrap_nrow = NULL,
                                ...) {
  x <- as.data.frame(x)
  # keep only columns needed for plotting
  if (!"ab" %in% colnames(x)) {
    x <- x[,c("gold_standard", "test", "essential_agreement")]
  } else {
    x <- x[,c("gold_standard", "test", "essential_agreement", "ab")]
  }

  if (match_axes) {
    x[["gold_standard"]] <- match_levels(x[["gold_standard"]], match_to = x[["test"]])
    x[["test"]] <- match_levels(x[["test"]], match_to = x[["gold_standard"]])
    if (add_missing_dilutions) {
      x[["gold_standard"]] <- fill_dilution_levels(x[["gold_standard"]])
      x[["test"]] <- fill_dilution_levels(x[["test"]])
    }
  }

  if (!"ab" %in% colnames(x) | length(unique(x[["ab"]])) == 1 |
      all(is.null(c(facet_wrap_ncol, facet_wrap_nrow)))) {
    p <- plot_mic_validation_single_ab(x, match_axes, ...)

    if ("ab" %in% colnames(x) & "mo" %in% colnames(x)) {
      bpoints <- AMR::clinical_breakpoints
      p <- p + ggplot2::geom_hline(yintercept = AMR::as.mic(bpoints[]))
    }
    p
  }
  else {
    p <- plot_mic_validation_multi_ab(x, match_axes = match_axes,
                                      facet_wrap_ncol = facet_wrap_ncol,
                                      facet_wrap_nrow = facet_wrap_nrow,
                                      ...)
    p
  }
}

#' Summary of MIC validation results
#'
#' @param object S3 mic_validation object
#' @param ... further optional parameters
#'
#' @export
#'
#' @return S3 mic_validation_summary object
#'
#' @description
#' Summarise the results of an MIC validation generated using compare_mic().
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' summary(val)
#' # or, for more detailed results
#' as.data.frame(summary(val))
summary.mic_validation <- function(object,
                                   ...) {
  if (!"ab" %in% names(object) & !"mo" %in% names(object)) {
    output <- list(EA_n = sum(object$essential_agreement == TRUE),
                   EA_pcent = sum(object$essential_agreement == TRUE) / length(object$essential_agreement),
                   bias = bias(object$gold_standard, object$test),
                   n = length(object$essential_agreement))
    class(output) <- append(class(output), "mic_validation_summary", after = 0)
    return(output)
  }

  if ("ab" %in% names(object) & !"mo" %in% names(object)) {
      output <- object |>
        as.data.frame() |>
        dplyr::group_by(.data[["ab"]]) |>
        dplyr::summarise(EA_n = sum(.data[["essential_agreement"]] == TRUE),
                         EA_pcent = sum(.data[["essential_agreement"]] == TRUE) / length(.data[["essential_agreement"]]),
                         bias = bias(object$gold_standard, object$test),
                         n = dplyr::n())
      class(output) <- append(class(output), "mic_validation_summary", after = 0)
      return(output)
  }

  if ("ab" %in% names(object) & "mo" %in% names(object)) {
    output <- object |>
        as.data.frame() |>
        dplyr::group_by(.data[["ab"]], .data[["mo"]]) |>
        dplyr::summarise(
          EA_pcent = sum(.data[["essential_agreement"]] == TRUE) / length(.data[["essential_agreement"]]),
          bias = bias(.data[["gold_standard"]], .data[["test"]]),
          resistant_pcent = AMR::proportion_R(.data[["gold_standard_sir"]], minimum = 1, as_percent = FALSE) * 100,
          minor_error_pcent = sum(.data[["error"]] == "m", na.rm = TRUE) / length(.data[["error"]]) * 100,
          major_error_pcent = sum(.data[["error"]] == "M", na.rm = TRUE) / length(.data[["error"]]) * 100,
          very_major_error_pcent = sum(.data[["error"]] == "vM", na.rm = TRUE) /length(.data[["error"]]) * 100,
          EA_n = sum(.data[["essential_agreement"]] == TRUE),
          resistant_n = AMR::count_resistant(.data[["gold_standard_sir"]]),
          minor_error_n = sum(.data[["error"]] == "m", na.rm = TRUE),
          major_error_n = sum(.data[["error"]] == "M", na.rm = TRUE),
          very_major_error_n = sum(.data[["error"]] == "vM", na.rm = TRUE),
          n = dplyr::n())
    class(output) <- append(class(output), "mic_validation_summary", after = 0)
    return(output)

  }
}

#' Print MIC validation summary
#'
#' @param x mic_validation_summary object
#' @param ... additional arguments
#'
#' @return character
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' print(summary(val))
print.mic_validation_summary <- function(x, ...) {
  if (!"ab" %in% names(x) & !"mo" %in% names(x)) {
    return(
      cat(glue::glue("MIC validation summary
                     Essential agreement: {x$EA_n} ({round(x$EA_pcent * 100, 2)}%)
                     Bias: {x$bias}"))
    )
  }
  if ("ab" %in% names(x) & !"mo" %in% names(x)) {
    antibiotic_str <- paste(unique(x$ab), collapse = ", ")
    return(
      cat(glue::glue("MIC validation summary
                     Number of observations: {sum(x$n)}
                     Antibiotic: {antibiotic_str}
                     Essential agreement: {sum(x$EA_n)} ({round(sum(x$EA_n) / sum(x$n) * 100, 2)}%)
                     Mean bias: {mean(x$bias)}
                     *Use as.data.frame() to see full summary*"))
    )
  }
  if (length(x) > 1 & "ab" %in% names(x) & "mo" %in% names(x)) {
    antibiotic_str <- paste(unique(x$ab), collapse = ", ")
    organism_str <- paste(unique(x$mo), collapse = ", ")

    return(
      cat(glue::glue("MIC validation summary
                     Antibiotic: {antibiotic_str}
                     Organism: {organism_str}
                     Essential agreement: {sum(x$EA_n)} ({round(sum(x$EA_n) / sum(x$n) * 100, 2)}%)
                     Resistant: {sum(x$resistant_n)} ({round(sum(x$resistant_n) / sum(x$n) * 100, 2)}%)
                     Minor errors: {sum(x$minor_error_n)} ({round(sum(x$minor_error_n) / sum(x$n) * 100, 2)}%)
                     Major errors: {sum(x$major_error_n)} ({round(sum(x$major_error_n) / sum(x$n) * 100, 2)}%)
                     Very major errors: {sum(x$very_major_error_n)} ({round(sum(x$very_major_error_n) / sum(x$n) * 100, 2)}%)
                     Mean bias: {mean(x$bias)}
                     N: {sum(x$n)}
                     *Use as.data.frame() to see full summary*"))
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
#'
#' @description
#' Check whether MIC values are within acceptable range for
#' quality control (QC). Every MIC experiment should include a control strain
#' with a known MIC. The results of the experiment are only valid if the control
#' strain MIC falls within the acceptable range. This function checks whether
#' an MIC result is within the acceptable range given: 1) a control strain
#' (usually identified as an ATCC or NCTC number), 2) an antibiotic name, and 3)
#' a guideline (EUCAST or CLSI). The acceptable range is defined by 'QC_table',
#' which is a dataset which is loaded with this package.
#'
#' The source of the QC values is the WHONET QC Ranges and Targets available from
#' the 'Antimicrobial Resistance Test Interpretation Engine' (AMRIE) repository:
#' https://github.com/AClark-WHONET/AMRIE
#'
#' @export
#' @references
#' O’Brien TF, Stelling JM. WHONET: An Information System for Monitoring
#' Antimicrobial Resistance. Emerg Infect Dis. 1995 Jun;1(2):66–66.

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
    qc_subset <- QC_table[QC_table$STRAIN == strain &
                            QC_table$ANTIBIOTIC == ab &
                            QC_table$GUIDELINE == guideline &
                            QC_table$YEAR == year &
                            QC_table$METHOD == "MIC", ]

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
#' `r lifecycle::badge('experimental')`
#'
#' @param measurement measured QC MIC
#' @param strain control strain identifier (usually ATCC)
#' @param ab antibiotic name (will be coerced to AMR::as.ab)
#' @param ignore_na ignores NA (returns TRUE)
#' @param guideline Guideline to use (EUCAST or CLSI)
#' @param year Guideline year (version)
#' @description
#' MIC experiments should include a control strain with a known MIC.
#' The MIC result for the control strain should be a particular target MIC. This
#' function checks whether the target MIC was achieved given: 1) a control
#' strain (usually identified as an ATCC or NCTC number), 2) an antibiotic name,
#' and 3) a guideline (EUCAST or CLSI).
#'
#' Since QC target values are currently not publicly available in an easy to
#' use format, this function takes a pragmatic approach -- for most antibiotics
#' and QC strains, the target is assumed to be the midpoint of the acceptable
#' range. This approximation is not necessarily equal to the QC target reported
#' by guideline setting bodies such as EUCAST. Therefore, this function is
#' considered experimental and should be used with caution.
#'
#' This function can be used alongnside qc_in_range(), which checks whether the
#' MIC is within the acceptable range.
#'
#' The source of the QC values is the WHONET QC Ranges and Targets available from
#' the 'Antimicrobial Resistance Test Interpretation Engine' (AMRIE) repository:
#' https://github.com/AClark-WHONET/AMRIE
#'
#' @return logical vector
#' @export
#'
#' @references
#' O’Brien TF, Stelling JM. WHONET: An Information System for Monitoring
#' Antimicrobial Resistance. Emerg Infect Dis. 1995 Jun;1(2):66–66.
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

    qc_subset <- QC_table[QC_table$STRAIN == strain &
                            QC_table$ANTIBIOTIC == ab &
                            QC_table$GUIDELINE == guideline &
                            QC_table$YEAR == year & QC_table$METHOD == "MIC", ]

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
#' `r lifecycle::badge('experimental')`
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
#' @description
#' MIC experiments are generally quality-controlled by including a control strain
#' with a known MIC. The MIC result for the control strain should be a particular
#' target MIC, or at least within an acceptable range. This function standardises
#' a measured MIC to the target MIC given: 1) a control strain (usually identified
#' as an ATCC or NCTC number), 2) an antibiotic name, and 3) a guideline (EUCAST
#' or CLSI). The definition of standardisation in this context is to adjust the
#' measured MIC based on the QC MIC. This is based on the following principles
#' and assumption:
#'
#' 1. A measured MIC is composed of two components: the true MIC and a
#' measurement error. The measurement error is considered to be inevitable when
#' measuring MICs, and is likely to be further composed of variability in
#' laboratory conditions and operator interpretation.
#' 2. It is assumed that the MIC of the control strain in the experiment has
#' also been affected by this error.
#'
#' The standardisation applied by this function uses the measured QC strain
#' MIC as a reference point, and scales the rest of the MICs to this reference.
#' In general, this means that the MICs are doubled or halved, depending on the
#' result of the QC MIC. A worked example is provided below and illustrates
#' the transformation that this function applies.
#'
#' There is no current evidence base for this approach, therefore, this
#' function is considered experimental and should be used with caution.
#'
#' @examples
#' # Ref strain QC MIC for GEN is 0.5
#' standardise_mic(
#'   test_measurement = c(AMR::as.mic(">8.0"),  # QC = 1, censored MIC remains censored
#'                        AMR::as.mic(4.0),  # QC = 0.5 which is on target, so stays same
#'                        AMR::as.mic(2),  # QC = 1, so scaled down to 1
#'                        AMR::as.mic(2)),  # QC = 0.25, so scaled up to 8
#'   qc_measurement = c(AMR::as.mic(1),
#'                      AMR::as.mic(0.5),
#'                      AMR::as.mic(1),
#'                      AMR::as.mic(0.25)),
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

#' Compare SIR results and generate categorical agreement
#'
#' @param gold_standard Susceptibility results in AMR::sir format
#' @param test Susceptibility results in AMR::sir format
#'
#' @return factor vector with the following levels: M, vM, m.
#' @export
#'
#' @description
#' Compare two AMR::sir vectors and generate a categorical agreement vector with
#' the following levels: M (major error), vM (very major error), m (minor error).
#' The error definitions are:
#'
#' 1. Major error (M): The test result is resistant (R) when the gold standard
#' is susceptible (S).
#' 2. vM (very major error): The test result is susceptible (S) when the gold
#' standard is resistant (R).
#' 3. Minor error (m): The test result is intermediate (I) when the gold standard
#' is susceptible (S) or resistant (R), or vice versa.
#'
#' @examples
#' gold_standard <- c("S", "R", "I", "I")
#' gold_standard <- AMR::as.sir(gold_standard)
#' test <- c("S", "I", "R", "R")
#' test <- AMR::as.sir(test)
#' compare_sir(gold_standard, test)
compare_sir <- function(gold_standard, test) {
  if (length(gold_standard) != length(test)) {
    stop("Inputs to compare_sir must be same length")
  }
  if (!all(AMR::is.sir(gold_standard)) | !all(AMR::is.sir(test))) {
    stop("Inputs to compare_sir must be AMR::sir objects")
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
      ),
      levels = c("M", "vM", "m")
    )
  )
}

#' Calculate MIC bias
#'
#' @param gold_standard AMR::mic vector
#' @param test AMR::mic vector
#'
#' @return numeric value
#' @export
#'
#' @description
#' Calculate the bias between two AMR::mic vectors. The bias is calculated as
#' the percentage of test MICs that are above the gold standard MICs minus the
#' percentage of test MICs that are below the gold standard MICs.
#'
#' @references
#' International Organization for Standardization. ISO 20776-2:2021
#' Available from: https://www.iso.org/standard/79377.html
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' bias(gold_standard, test)
bias <- function(gold_standard, test) {
  if (length(gold_standard) != length(test)) {
    stop("Inputs to bias must be same length")
  }
  bias_above <- sum(test > gold_standard)
  bias_below <- sum(test < gold_standard)

  n <- length(gold_standard)
  return(bias_above / n * 100 - bias_below / n * 100)
}

#' Table
#'
#' @param x an mic_validation object to be tabulated into an essential agreement
#' frequency table, or object/s to be passed to base::table
#' @param ... further arguments
#'
#' @rdname table
#' @export
#'
table <- function(x, ...) {
  UseMethod("table")
}

#' @rdname table
#' @export
table.default <- function(x, ...) {
  return(base::table(x, ...))
}

tabulate_flex <- function(t, ea, bold, ea_color, gold_standard_name, test_name) {
  t_flex <- t |>
    stats::addmargins(FUN = list(Total = sum),
                      quiet = TRUE) |>
    as.data.frame.matrix() |>
    tibble::rownames_to_column(var = "test_mics") |>
    dplyr::mutate("test" = test_name, .before = "test_mics") |>
    flextable::flextable()

  for (i in 1:nrow(t)) {
    for (j in 1:ncol(t)) {
      if (ea[i, j] == TRUE) {
        if (bold) {
          t_flex <- flextable::bold(t_flex, i = i, j = j + 2)
        }
        if (!is.null(ea_color)) {
          t_flex <- flextable::bg(t_flex, i = i, j = j + 2, bg = ea_color)
        }
      }
    }
  }
  t_flex <- t_flex |>
    flextable::set_header_labels("test_mics" = "", test = "") |>
    flextable::add_header_row(values = c("", gold_standard_name, ""),
                              colwidths = c(1, length(t_flex[["col_keys"]]) - 2, 1)) |>
    flextable::align(align = "center", i=1,j=2,part = "header") |>
    flextable::align(align = "right", j=2) |>
    flextable::merge_v(j = 1) |>
    flextable::border_remove() |>
    flextable::align(align = "right", part = "body")
  t_flex
}

#' @rdname table
#'
#' @param x mic_validation S3 object
#' @param format simple or flextable
#' @param fill_dilutions Fill dilutions that are not present in the data in
#' order to match the y- and x- axes
#' @param bold Bold cells where essential agreement is TRUE
#' @param ea_color Background color for essential agreement cells
#' @param gold_standard_name Name of the gold standard to display in output
#' @param test_name Name of the test to display in output
#' @param ... further arguments
#'
#' @return table or flextable object
#'
#' @export
#'
#' @examples
#' gold_standard <- c("<0.25", "8", "64", ">64")
#' test <- c("<0.25", "2", "16", "64")
#' val <- compare_mic(gold_standard, test)
#' table(val)
table.mic_validation <- function(x,
                                 format = 'flextable',
                                 fill_dilutions = TRUE,
                                 bold = TRUE,
                                 ea_color = NULL,
                                 gold_standard_name = "Gold Standard",
                                 test_name = "Test",
                                 ...) {
  test <- match_levels(x[["test"]], match_to = x[["gold_standard"]])
  gold_standard <- match_levels(x[["gold_standard"]], match_to = x[["test"]])

  if (fill_dilutions) {
    test <- fill_dilution_levels(test)
    gold_standard <- fill_dilution_levels(gold_standard)
  }

  t <- table(test,
             gold_standard)

  rnames <- rownames(t)
  cnames <- colnames(t)
  grid <- expand.grid(rnames, cnames)
  ea <- essential_agreement(as.character(grid$Var1), as.character(grid$Var2))
  ea <- matrix(ea, nrow = length(rnames))
  ea <- as.data.frame(ea)

  if (format == "simple") {
    return(t)
  }
  if (format == "flextable") {
    return(tabulate_flex(t, ea, bold, ea_color, gold_standard_name, test_name))
  }

  stop("Format must be simple or flextable")
}
