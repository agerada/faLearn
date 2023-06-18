kmer_data_to_csv <- function(data, target_path, overwrite = FALSE) {
  # Data must be in the following format:
  # dataframe (or equivalent) with cols
  # genome_id, V1 .. Vn (for n kmers)

  if (file.exists(target_path) && !overwrite) {
    message("Skipping backup as file is already present. 
    Force overwrite if required")
    return()
  }

  if (file.exists(target_path)) {
    message(paste("Overwriting kmer backup to", target_path))
    write.csv(data, target_path, row.names = FALSE)
    return()
  }

  message(paste("Saving kmer backup to", target_path))
  write.csv(data, target_path, row.names = FALSE)
}
