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

usethis::use_data(QC_table, overwrite = TRUE, internal = TRUE)
