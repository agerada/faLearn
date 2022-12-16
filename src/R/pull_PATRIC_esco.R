library(tidyverse)
library(parallel)
library(glue)
library(Biostrings)
library(Rcpp)
library(bench)

# directories
database_txt_dir <- "data/databases/patric/PATRIC_genomes_AMR.txt"
genomes_download_dir <- "data/genomes/patric/"

# make list of unique genome_id for download
patric_amr_list <- read_delim(database_txt_dir, delim='\t')

# select only samples with broth or agar dilution
patric_amr_list_esco_agar_broth <- patric_amr_list %>% 
  filter(str_detect(genome_name, "Escherichia coli")) %>% 
  filter(measurement_unit == "mg/L")

esco_genome_ids <- unique(patric_amr_list_esco_agar_broth$genome_id)

# set download dir
old_wd <- getwd()
setwd(genomes_download_dir)

n_downloads <- length(esco_genome_ids)
# check we are not trying to downloading more than available
n_downloads <- ifelse(length(esco_genome_ids) < n_downloads, 
                      length(esco_genome_ids), n_downloads)

esco_genome_paths <- glue("ftp://ftp.patricbrc.org/genomes/{esco_genome_ids}/{esco_genome_ids}.fna")

i <- 1
while(i <= n_downloads){
  print(glue("Downloading file {i} of {n_downloads}"))
  tryCatch(download.file(esco_genome_paths[[i]], 
                destfile = glue("{esco_genome_ids[[i]]}.fna"), 
                mode="wb"), 
           error = function(e) print(glue("Unable to download {esco_genome_ids[[i]]}"))
           )
  i <- i+1
}
setwd(old_wd)
