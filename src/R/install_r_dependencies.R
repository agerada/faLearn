#!/usr/bin/env Rscript
# Copyright 2022 Alessandro Gerada alessandro.gerada@liverpool.ac.uk

packages <- c("tidyverse", "parallel", "glue", "Rcpp", "snow", "pbapply")
installed <- installed.packages()

if(!("BiocManager" %in% installed)) install.packages("BiocManager")
# install core Bioconductor
message("Installing/updating Bioconductor")
BiocManager::install()
# install Biostrings
message("Installing Biostrings")
BiocManager::install("Biostrings")

install.packages(setdiff(packages, installed))

# test Rcpp
message("Testing Rcpp")
if (Rcpp::evalCpp("2+2") == 4) {
  message("Rcpp functioning correctly.")
} else {
  stop("Rcpp not functioning. Check installation and C++ compiler. ")
}
