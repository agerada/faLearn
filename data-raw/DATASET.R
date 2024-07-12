## code to prepare QC_table dataset goes here

QC_table <- read.csv("data-raw/eucast_qc.csv")
QC_table$ab <- AMR::as.ab(QC_table$ab)
QC_table <- QC_table %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("mic_"), AMR::as.mic))

usethis::use_data(QC_table, overwrite = TRUE, internal = TRUE)
