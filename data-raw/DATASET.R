library(devtools)
devtools::load_all()

## Code to prepare qc_targets dataset goes here

## The source of this dataset is the WHONET QC Ranges and Targets available from
## the 'Antimicrobial Resistance Test Interpretation Engine' (AMRIE) repository:
## https://github.com/AClark-WHONET/AMRIE

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
n <- 100
ab <- AMR::as.ab("Gentamicin")
mo <- AMR::as.mo("Escherichia coli")
bimodal_probs <- c(1,1,2,3,4,5,4,3,2,1,1,
                   1,2,3,4,3,2,1,1)
# bimodal_probs <- bimodal_probs / sum(bimodal_probs)
gs <- sample(mic_range(), size = n, prob = bimodal_probs, replace = T)
gs <- AMR::as.mic(gs)

# now to generate the test values, assume that the value is the same as gs
# but with added noise
# log2 transform needed when adding noise otherwise lower MICs get proportionately
# more noise
test_values <- log2(gs) + rnorm(length(gs), 0, 1)
test_values <- 2^test_values
test_values[test_values < min(mic_range())] <- min(mic_range())
test_values <- force_mic(test_values)
test_values <- AMR::as.mic(test_values)

example_mics <- data.frame(gs = gs,
                           test_values = test_values)
example_mics$mo <- mo
example_mics$ab <- ab

readr::write_tsv(example_mics, "data-raw/example_mics.txt")

use_data(QC_table, ecoffs, example_mics,
         overwrite = TRUE, internal = TRUE)
