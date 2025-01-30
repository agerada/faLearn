library(devtools)
devtools::load_all()

## Code to prepare qc_targets dataset goes here

## The source of this dataset is the WHONET QC Ranges and Targets available from
## the 'Antimicrobial Resistance Test Interpretation Engine' (AMRIE) repository:
## https://github.com/AClark-WHONET/AMRIE which is made available under the
## GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007.

whonet_qc_ranges <- readr::read_tsv(
  "https://raw.githubusercontent.com/AClark-WHONET/AMRIE/main/Interpretation%20Engine/Resources/QC_Ranges.txt",
  col_types = readr::cols(.default = "c"))
readr::write_tsv(whonet_qc_ranges, "data-raw/whonet_qc_ranges.txt")

qc_targets <- readr::read_csv("data-raw/qc_targets.csv",
                              col_types = readr::cols(.default = "c"))
QC_table <- merge(whonet_qc_ranges,
                  qc_targets,
                  by = c("GUIDELINE", "YEAR", "STRAIN", "WHONET_ABX_CODE", "METHOD"),
                  all.x = TRUE)

QC_table$ANTIBIOTIC <- AMR::as.ab(QC_table$ANTIBIOTIC)

## where min and max target are missing from qc_targets.csv, use the mid
## point MIC between minimum and maximum instead

mid_mic <- function(lower, upper) {
  diff <- AMR::as.mic(upper) - AMR::as.mic(lower)
  (diff / 2) + AMR::as.mic(lower)
}

QC_table$TARGET_ESTIMATED <- ifelse(is.na(QC_table$MINIMUM_TARGET) | is.na(QC_table$MAXIMUM_TARGET),
                                    TRUE,
                                    FALSE)

QC_table$MINIMUM_TARGET <- ifelse(
  is.na(QC_table$MINIMUM_TARGET) & QC_table$METHOD == "MIC",
  mid_mic(QC_table$MINIMUM, QC_table$MAXIMUM),
  QC_table$MINIMUM_TARGET)

QC_table$MAXIMUM_TARGET <- ifelse(
  is.na(QC_table$MAXIMUM_TARGET) & QC_table$METHOD == "MIC",
  QC_table$MINIMUM_TARGET,
  QC_table$MAXIMUM_TARGET)

## repeat for disks

mid_disk <- function(lower, upper) {
  diff <- AMR::as.disk(upper) - AMR::as.disk(lower)
  (diff / 2 ) + AMR::as.disk(lower)
}

QC_table$MINIMUM_TARGET <- ifelse(
  is.na(QC_table$MINIMUM_TARGET) & QC_table$METHOD == "Disk",
  mid_disk(QC_table$MINIMUM, QC_table$MAXIMUM),
  QC_table$MINIMUM_TARGET)

QC_table$MAXIMUM_TARGET <- ifelse(
  is.na(QC_table$MAXIMUM_TARGET) & QC_table$METHOD == "Disk",
  QC_table$MINIMUM_TARGET,
  QC_table$MAXIMUM_TARGET)


## ECOFFS

ecoffs_files <- list.files("data-raw/ecoffs", full.names = TRUE)
ecoffs <- readr::read_csv(ecoffs_files,
                          col_types = readr::cols(.default = "c"),
                          id = "organism") |>
  dplyr::mutate(organism = basename(organism)) |>
  dplyr::mutate(organism = tools::file_path_sans_ext(organism))
ecoffs <- ecoffs |>
  dplyr::rename(antibiotic = ...1) |>
  dplyr::mutate(antibiotic = AMR::as.ab(antibiotic))

ecoffs$organism <- AMR::as.mo(ecoffs$organism)

## Example dataset
set.seed(42)
n <- 300
ab <- c(rep(AMR::as.ab("Gentamicin"), n / 3),
        rep(AMR::as.ab("Meropenem"), n / 3),
        rep(AMR::as.ab("Amoxicillin"), n / 3))
mo <- rep(AMR::as.mo("Escherichia coli"), n)
bimodal_probs <- c(1,1,2,3,4,5,4,3,2,1,1,
                   1,2,3,4,3,2,1,1)
bimodal_probs <- bimodal_probs / sum(bimodal_probs)
res_probs <- c(0.01,0.01,2,3,4,5,4,3,2,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
res_probs <- res_probs / sum(res_probs)
sens_probs <- c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,2,3,4,3,2,0.01,0.01)
sens_probs <- sens_probs / sum(sens_probs)

gs <- sample(mic_range(), size = n / 3, prob = bimodal_probs, replace = T)
gs <- AMR::as.mic(gs)
# now to generate the test values, assume that the value is the same as gs
# but with added noise
# log2 transform needed when adding noise otherwise lower MICs get proportionately
# more noise
bimodal_values <- log2(gs) + rnorm(length(gs), 0, 1)
bimodal_values <- 2^bimodal_values
bimodal_values[bimodal_values < min(mic_range())] <- min(mic_range())
bimodal_values <- force_mic(bimodal_values)
bimodal_values <- AMR::as.mic(bimodal_values)

# repeat for sensitivity and resistance
sens <- sample(mic_range(), size = n / 3, prob = sens_probs, replace = T)
sens <- AMR::as.mic(sens)
sens_values <- log2(sens) + rnorm(length(sens), 0, 1)
sens_values <- 2^sens_values
sens_values[sens_values < min(mic_range())] <- min(mic_range())
sens_values <- force_mic(sens_values)
sens_values <- AMR::as.mic(sens_values)

res <- sample(mic_range(), size = n / 3, prob = res_probs, replace = T)
res <- AMR::as.mic(res)
res_values <- log2(res) + rnorm(length(res), 0, 1)
res_values <- 2^res_values
res_values[res_values < min(mic_range())] <- min(mic_range())
res_values <- force_mic(res_values)
res_values <- AMR::as.mic(res_values)

# combine the values
values <- c(bimodal_values, sens_values, res_values)
gs <- c(gs, sens, res)

example_mics <- data.frame(gs = gs, test = values)
example_mics$mo <- mo
example_mics$ab <- ab

readr::write_tsv(example_mics, "data-raw/example_mics.txt")

use_data(QC_table, overwrite = TRUE, internal = TRUE)
use_data(example_mics, ecoffs, overwrite = TRUE, internal = FALSE)
