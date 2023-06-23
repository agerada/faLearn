#!/usr/bin/env Rscript
# Copyright 2022 Alessandro Gerada alessandro.gerada@liverpool.ac.uk

packages <- c("tidyverse", "parallel", "glue", "Rcpp", "snow", "pbapply",
              "optparse", "here",
              "AMR", "multidplyr")
installed <- installed.packages()

if(!("BiocManager" %in% installed)) install.packages("BiocManager", repos = "http://cran.us.r-project.org")
# install core Bioconductor
message("Installing/updating Bioconductor")
BiocManager::install()
# install Biostrings
message("Installing Biostrings")
BiocManager::install("Biostrings")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ShortRead")

#
install.packages(setdiff(packages, installed), repos = "http://cran.us.r-project.org")

# test Rcpp
message("Testing Rcpp")
if (Rcpp::evalCpp("2+2") == 4) {
  message("Rcpp functioning correctly.")
} else {
  stop("Rcpp not functioning. Check installation and C++ compiler. ")
}
