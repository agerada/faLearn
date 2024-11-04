## code to prepare qc_targets dataset goes here

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

ecoffs <- readr::read_csv("data-raw/ecoffs.csv",
                          col_types = readr::cols(.default = "c"))
ecoffs <- ecoffs %>%
  dplyr::rename(antibiotic = ...1) %>%
  dplyr::mutate(antibiotic = AMR::as.ab(antibiotic))

## Example dataset

example_raw <- readr::read_tsv("data-raw/example_raw.txt",
                                col_types = readr::cols(.default = "c")) %>%
  dplyr::mutate(genome_name = AMR::mo_name(AMR_org)) %>%
  dplyr::relocate(genome_id,
                  date_collected,
                  AMR_org,
                  antibiotic,
                  measurement,
                  measurement_unit,
                  laboratory_typing_method,
                  resistant_phenotype,
                  lca_class,
                  dplyr::starts_with("POS_QC"),
                  dplyr::starts_with("QC"))

# get the QC MIC (E. coli ATC 25922) for each antimicrobial
qc_mic <- example_raw %>%
  dplyr::distinct(genome_id, .keep_all = TRUE) %>%
  dplyr::select(genome_id, dplyr::starts_with("POS_QC_")) %>%
  tidyr::pivot_longer(dplyr::starts_with("POS_QC_"), names_to = "antibiotic",
               values_to = "qc_mic") %>%
  dplyr::mutate(antibiotic = stringr::str_remove(antibiotic, "POS_QC_")) %>%
  dplyr::mutate(laboratory_typing_method = "Agar dilution")

# get the growth control for each antimicrobial (P/W/F)
qc_control <- example_raw %>%
  dplyr::distinct(genome_id, .keep_all = TRUE) %>%
  dplyr::select(genome_id, dplyr::starts_with("QC_")) %>%
  tidyr::pivot_longer(dplyr::starts_with("QC_"), names_to = "antibiotic",
               values_to = "qc_growth") %>%
  dplyr::mutate(antibiotic = stringr::str_remove(antibiotic, "QC_")) %>%
  dplyr::mutate(laboratory_typing_method = "Agar dilution")

# Combine and remove any antimicrobials with QC outside range
example_mics <- example_raw %>%
  dplyr::left_join(qc_mic, by = c("genome_id", "antibiotic", "laboratory_typing_method")) %>%
  dplyr::left_join(qc_control, by = c("genome_id", "antibiotic", "laboratory_typing_method")) %>%
  dplyr::filter(laboratory_typing_method == "Agar dilution") %>%
  dplyr::mutate(measurement = AMR::as.mic(measurement)) %>%
  dplyr::mutate(qc_mic = AMR::as.mic(qc_mic)) %>%
  dplyr::filter(qc_in_range(qc_mic, 25922, antibiotic)) %>%
  dplyr::filter(qc_growth != "F") %>%
  dplyr::select(!dplyr::starts_with("QC_") & !dplyr::starts_with("POS_QC_"))

example_mics <- example_mics %>%
  as_patric_db()

example_discs <- example_raw %>%
  dplyr::filter(laboratory_typing_method == "Disk diffusion") %>%
  dplyr::mutate(measurement = AMR::as.disk(measurement)) %>%
  dplyr::select(!dplyr::starts_with("QC_") & !dplyr::starts_with("POS_QC_"))

example_discs <- example_discs %>%
  as_patric_db()

usethis::use_data(QC_table, ecoffs,
  example_mics, example_discs, overwrite = TRUE, internal = TRUE)
