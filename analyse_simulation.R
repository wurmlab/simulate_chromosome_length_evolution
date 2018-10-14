#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-i", "--insertion_rate"), type="numeric", default=5e-6,
              help="Insertion mutation rate [default %default]",
              metavar="number"),
  make_option(c("-d", "--deletion_rate"), type="numeric", default=1e-5,
              help="Insertion mutation rate [default %default]",
              metavar="number"),
  make_option(c("-I", "--intactness_matrix"), type="character", default="int_out.txt",
              help="File for output intactness matrix [default= %default]", metavar="character"),
  make_option(c("-c", "--chromosome_list_dir"), type="character", default="chr_out",
              help="Directory for output chromosome list [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)


# Run setup
generation_number  <- 1000
# Chromosome setup
N                  <- 1000 # number of individuals
gene_length        <- 2000
intergene_length   <- 6000
gene_number        <- 100

# Output
intactness_matrix_output <- opt$intactness_matrix
chromosome_list_output   <- opt$chromosome_list_dir

# Save which
every_nth <- 20

# ---------------------------------------------------------------------------- #
# Saved generations
save_every_nth_seq <- seq(from = every_nth, to = N, by = every_nth)


# ---------------------------------------------------------------------------- #
# Open matrix

intactness_matrix <- read.table(intactness_matrix_output)

# Divide into each generation
generation_boundaries <- seq(from = 1, to = N * length(save_every_nth_seq), by = N)
intactness_per_gen    <- lapply(generation_boundaries,
                                function(i) intactness_matrix[i:(i + N - 1), ])

# ---------------------------------------------------------------------------- #
# Evolution of intactness per generation

generationChromosomeMeanIntactness <- function(intactness_matrix) {
  mean_intactness_per_chromosome <- apply(intactness_matrix, 1, mean)
}

mean_intactness_perchr_pergen <- lapply(intactness_per_gen,
                                        generationChromosomeMeanIntactness)

mean_mean_intactness_perchr_pergen <- sapply(mean_intactness_perchr_pergen, mean)

numberOfBrokenGenes <- function(intactness_matrix) {
  return(apply(intactness_matrix, 1, function(v) sum(v == 0)))
}

broken_genes_pergen <- lapply(intactness_per_gen,
                              numberOfBrokenGenes)
mean_broken_genes_pergen <- sapply(broken_genes_pergen, mean)

# ---------------------------------------------------------------------------- #
# Evolution size per generation

chrom_list_out <- paste(chromosome_list_output, save_every_nth_seq, sep="/")
chrom_list_out <- paste(chrom_list_out, ".Rda", sep= "")

chr_lengths       <- list()
gene_lengths      <- list()
intergene_lengths <- list()

chr_rle_list <- list()

for (i in 1:length(chrom_list_out)) {
  print(paste("Loading from", chrom_list_out[[i]]))
  load(chrom_list_out[[i]])

  chr_lengths[[i]] <- sapply(chromosome_list, length)
  gene_lengths[[i]]   <- sapply(chromosome_list, function(y) sum(y != 0))
  intergene_lengths[[i]] <- sapply(chromosome_list, function(y) sum(y == 0))

  # chr_rle_list[[i]] <- lapply(chromosome_list, rle)
}

chr_lengths <- lapply(chr_list, function(x) sapply(x, length))

max_chr_length    <- lapply(chr_lengths, max)
min_chr_length    <- lapply(chr_lengths, min)
mean_chr_length   <- lapply(chr_lengths, mean)
sd_chr_length     <- lapply(chr_lengths, sd)
median_chr_length <- lapply(chr_lengths, median)
mad_chr_length    <- lapply(chr_lengths, mad)

gene_lengths <- lapply(chr_list, function(x) sapply(x, function(y) sum(y != 0)))
max_gene_length    <- lapply(gene_lengths, max)
min_gene_length    <- lapply(gene_lengths, min)
mean_gene_length   <- lapply(gene_lengths, mean)
sd_gene_length     <- lapply(gene_lengths, sd)
median_gene_length <- lapply(gene_lengths, median)
mad_gene_length    <- lapply(gene_lengths, mad)

intergene_lengths <- lapply(chr_list, function(x) sapply(x, function(y) sum(y == 0)))
max_intergene_length    <- lapply(intergene_lengths, max)
min_intergene_length    <- lapply(intergene_lengths, min)
mean_intergene_length   <- lapply(intergene_lengths, mean)
sd_intergene_length     <- lapply(intergene_lengths, sd)
median_intergene_length <- lapply(intergene_lengths, median)
mad_intergene_length    <- lapply(intergene_lengths, mad)
