backup_kmer_data <- function(
  data, target_path, overwrite = FALSE,
  file_format = "csv") {
  # Data must be in the following format:
  # dataframe (or equivalent) with cols
  # genome_id, V1 .. Vn (for n kmers)
  # supports csv or parquet format

  writing_function <- switch(
    file_format,
    csv = function(x, path) {
      write.csv(x, path, row.names = FALSE)
    },
    parquet = function(x, path) {
      arrow::write_parquet(x, path)
    }
  )
  if (file.exists(target_path) && !overwrite) {
    message("Skipping backup as file is already present. 
    Force overwrite if required")
    return()
  }

  if (file.exists(target_path)) {
    message(paste("Overwriting kmer backup to", target_path))
    writing_function(data, target_path)
    return()
  }

  message(paste("Saving kmer backup to", target_path))
  writing_function(data, target_path)
}
