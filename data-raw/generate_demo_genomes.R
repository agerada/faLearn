#!/usr/bin/env Rscript
# Script to generate a demonstration dataset of genomes and MIC labels
# This script is reproducible (fixed seed) and will save the dataset to
# data/demo_genomes.rda when run from the package root.

## Packages
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Biostrings is required to run this script. Install with BiocManager::install('Biostrings')")
}

set.seed(42)

library(Biostrings)

# Parameters
n_genomes <- 100L
target_genome_size <- 5000L  # approx total bases per genome
min_contigs <- 1L
max_contigs <- 8L

# define two predictive gene sequences (short ~500bp genes)
make_random_gene <- function(len = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  paste0(sample(c("A","C","G","T"), len, replace = TRUE), collapse = "")
}

gene1_seq <- make_random_gene(20, seed = 123)
gene2_seq <- make_random_gene(25, seed = 456)

# choose genomes that receive each gene
gene1_idx <- sample(seq_len(n_genomes), size = 0.15 * n_genomes)
gene2_idx <- sample(seq_len(n_genomes), size = 0.10 * n_genomes)

## Ensure some overlap
both_idx <- sample(gene1_idx, size = 0.025 * n_genomes)
gene2_idx <- unique(c(gene2_idx, both_idx))

genomes <- vector("list", n_genomes)
names(genomes) <- paste0("ECOLI_", sprintf("%03d", seq_len(n_genomes)))

for (i in seq_len(n_genomes)) {
  n_contigs <- sample(min_contigs:max_contigs, 1)
  # random partition of target_genome_size across contigs
  cuts <- sort(sample(seq_len(target_genome_size - 1), n_contigs - 1))
  lens <- diff(c(0L, cuts, target_genome_size))

  contigs <- character(length(lens))
  for (j in seq_along(lens)) {
    contigs[j] <- paste0(sample(c("A","C","G","T"), lens[j], replace = TRUE), collapse = "")
  }

  # inject gene1 or gene2 into a random contig for selected genomes
  if (i %in% gene1_idx) {
    c_idx <- sample(seq_along(contigs), 1)
    insert_pos <- sample(seq_len(nchar(contigs[c_idx]) + 1), 1)
    contigs[c_idx] <- paste0(substr(contigs[c_idx], 1, insert_pos - 1),
                             gene1_seq,
                             substr(contigs[c_idx], insert_pos, nchar(contigs[c_idx])))
  }
  if (i %in% gene2_idx) {
    c_idx <- sample(seq_along(contigs), 1)
    insert_pos <- sample(seq_len(nchar(contigs[c_idx]) + 1), 1)
    contigs[c_idx] <- paste0(substr(contigs[c_idx], 1, insert_pos - 1),
                             gene2_seq,
                             substr(contigs[c_idx], insert_pos, nchar(contigs[c_idx])))
  }

  ds <- DNAStringSet(contigs)
  names(ds) <- paste0(names(genomes)[i], "_contig", seq_along(ds))
  genomes[[i]] <- ds
}

# create presence matrix (no genome names): two columns (gene1, gene2)
gene_presence_df <- data.frame(
  gene1 = names(genomes) %in% names(genomes)[gene1_idx],
  gene2 = names(genomes) %in% names(genomes)[gene2_idx],
  stringsAsFactors = FALSE
)
gene_presence <- as.matrix(gene_presence_df)
colnames(gene_presence) <- c("gene1", "gene2")
rownames(gene_presence) <- NULL

# generate MIC labels (mg/L). Wild-type: ~0.5, gene1 raises to ~4, gene2 to ~8,
# both to ~16; small noise added.
base_mic <- 0.5
mic <- numeric(n_genomes)
for (i in seq_len(n_genomes)) {
  g1 <- gene_presence[i, "gene1"]
  g2 <- gene_presence[i, "gene2"]
  if (!g1 && !g2) {
    m <- base_mic + rnorm(1, sd = 0.1)
  } else if (g1 && !g2) {
    m <- 4 + rnorm(1, sd = 0.2)
  } else if (!g1 && g2) {
    m <- 32 + rnorm(1, sd = 0.5)
  } else {
    m <- 128 + rnorm(1, sd = 1)
  }
  mic[i] <- round(pmax(m, 0.03), digits = 3)
}

## Package-friendly single object `example_genomes`
example_genomes <- list(
  genomes = genomes,                 # named list of Biostrings::DNAStringSet
  annotations = gene_presence,       # 2-column logical matrix (gene1, gene2) matching genomes order
  labels = log2(mic)                 # MIC labels on log2 scale (log2(mg/L))
)

# Save to data/ (when run from package root)
if (!dir.exists("data")) dir.create("data")
save(example_genomes, file = "data/example_genomes.rda")

message("Saved data/example_genomes.rda containing object: example_genomes (genomes, annotations, labels)")
