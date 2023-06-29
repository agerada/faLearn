library(Biostrings)
library(arrow)
library(data.table)

kmer_width <- 12
n_process <- 5
path <- "/volatile/agerada/molecularMIC/genomes/e_coli_all"
out_path <- "/volatile/agerada/molecularMIC/kmers/e_coli_all/arrow/"

fasta_paths <- list.files(path, full.names = TRUE)

fasta_paths <- head(fasta_paths, n = n_process)

cl <- parallel::makeCluster(56)
parallel::clusterExport(cl, c("fasta_paths", "oligonucleotideFrequency", "kmer_width"))

genomes <- pbapply::pblapply(fasta_paths, readDNAStringSet, cl = cl)

bench::mark(
    kmers <- pbapply::pblapply(
        genomes,
        function(x) {
            oligonucleotideFrequency(
                unlist(x), width = kmer_width, with.labels = FALSE)
            },
        cl = cl)
)

bench::mark({
    setDT(kmers)

}
)

bench::mark(
    kmers <- data.table::transpose(kmers)
)

arrow::write_parquet(kmers, file.path(out_path, "test.parquet"))
