read_kmer_r <- function(x, kmer = 3){ 
  n <- nchar(x) - kmer + 2
  i <- 1
  kmer_dict <- new.env(hash=TRUE)
  while (i < n){
    kmer_i <- substr(x, i, i + kmer - 1)
    if (kmer_i %in% names(kmer_dict)) {
      kmer_dict[[kmer_i]] <- kmer_dict[[kmer_i]] + 1
    } else {
      kmer_dict[[kmer_i]] <- 1
    }
    i <- i+1
  }
  return(as.list(kmer_dict))
}
