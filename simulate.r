#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------- #
# DESCRIPTION HERE



# ---------------------------------------------------------------------------- #

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


# ---------------------------------------------------------------------------- #
# DEFAULT SETTINGS

# Run setup
generation_number  <- 1000
# Chromosome setup
N                  <- 1000 # number of individuals
gene_length        <- 2000
intergene_length   <- 6000
gene_number        <- 100

# Mutation setup
# per base pair per individual per generation
insertion_rate     <- opt$insertion_rate
deletion_rate      <- opt$deletion_rate
mean_mutation_size <- 1000

# Cost of mutations
cost_of_intergenic_insertion <- 0.1
cost_of_intergenic_deletion  <- 0.1
minimum_intactness           <- 0.3

# Output
intactness_matrix_output <- opt$intactness_matrix
chromosome_list_output   <- opt$chromosome_list_dir

# Save which
every_nth <- 20

# ----------------------------------------------------------------------------------------- #

dir.create(chromosome_list_output, recursive=TRUE)

if (file.exists(intactness_matrix_output)) {
  stop(paste("File", intactness_matrix_output, "already exists. Stopping"))
}

save_every_nth_seq <- seq(from = every_nth, to = N, by = every_nth)

# ----------------------------------------------------------------------------------------- #

# Functions used in simulation

deletionInWhichGenes <- function(chromosome_coords, del_start, del_end) {
  deleted_region <- chromosome_coords[del_start:del_end]
  return(unique(deleted_region[deleted_region != 0]))
}

intergenicMutationEffect <- function(chromosome_coords, intactness, cost_of_intergenic_mutation, mutation_coord) {
  # mutation_coord needs to be the insertion start for insertions and deletion end for deletions
  # If mutation_coord falls in intergenic region (where genotype == 0)
  if (chromosome_coords[mutation_coord] == 0) {

    # Find next gene: the first value in coord vector after mutation_coord that is not 0
       # with an exception if there are no intergenic regions after the deletion
    possible_chromosome <- chromosome_coords[mutation_coord:length(chromosome_coords)]
    gene_positions      <- possible_chromosome != 0

    if (any(gene_positions)) {

      next_gene_along <- min(possible_chromosome[gene_positions])

      # Apply cost of mutation to intactness value
      new_intactness  <- intactness[next_gene_along] - cost_of_intergenic_mutation

      # Make sure intactness does not go below zero
      new_intactness[new_intactness < 0] <- 0

      # Update intacness vecor
      intactness[next_gene_along] <- new_intactness
    }
  }
  return(intactness)
}

mutationParamters <- function(chromosome_length, mutation_rate) {

  ## generate the number of mutations
  mutation_number   <- rpois(n = 1,
                             lambda = mutation_rate * chromosome_length)

  ## Mutation positions
  mutation_position <- sample(x       = 1:chromosome_length,
                              size    = mutation_number,
                              replace = FALSE)
  mutation_position <- sort(mutation_position)

  ## Mutation size
  mutation_size <- rpois(n      = mutation_number,
                         lambda = mean_mutation_size)
  return(list(position = mutation_position, size = mutation_size))

}


simulateDeletions <- function(chromosome,
                              intactness_vector,
                              cost_of_intergenic_mutation,
                              mutation_rate,
                              minimum_intactness) {

  # Set up the mutations
  deletions <- mutationParamters(chromosome_length = length(chromosome),
                                 mutation_rate = mutation_rate)
  deletions_position <- deletions$position
  deletions_size     <- deletions$size
  deletions_number   <- length(deletions_position)
  deletions_end      <- deletions_position + deletions_size

  # Edge case
  deletions_end[deletions_end > length(chromosome)] <- length(chromosome)

  if (deletions_number == 0) {
    survival = TRUE
  } else if (deletions_number > 0) {
    # Does any deletion affect a gene?
    deletions_in_genes <- lapply(1:deletions_number,
      function(i) deletionInWhichGenes(chromosome_coords = chromosome,
                                       del_start         = deletions_position[i],
                                       del_end           = deletions_end[i]))
    genes_deleted <- unique(do.call(c, deletions_in_genes))

    # only survive if all genes with deletion are no longer intact
    survival <- all(intactness_vector[genes_deleted] < minimum_intactness)

    # Don't chromosome if the individual dies anyway
    if (survival == TRUE) {
      # Deleted genes need to get no intactness left
      intactness_vector[genes_deleted] <- -9

      for (i in 1:deletions_number) {
        # Check whether there are deletions in intergenic regions and update the `intactness` vector
        intactness_vector <- intergenicMutationEffect(chromosome_coords = chromosome,
                    intactness                  = intactness_vector,
                    cost_of_intergenic_mutation = cost_of_intergenic_mutation,
                    mutation_coord              = deletions_end[i])
      }

      # Delete the positions from the chromosome
      for (i in 1:deletions_number) {
        chromosome <- chromosome[-c(deletions_position[i]:deletions_end[i])]
      }
    }
  } else {
    stop("Negative number of deletions??")
  }
   return(list(survival          = survival,
               chromosome        = chromosome,
               intactness_vector = intactness_vector))
}

simulateInsertions <- function(chromosome,
                               intactness_vector,
                               cost_of_intergenic_mutation,
                               mutation_rate,
                               minimum_intactness) {
  # Set up the mutations
  insertions          <- mutationParamters(chromosome_length = length(chromosome),
                                           mutation_rate     = mutation_rate)
  insertions_position <- insertions$position
  insertions_size     <- insertions$size
  insertions_number   <- length(insertions_position)

  if (insertions_number == 0) {
    survival = TRUE
  } else if (insertions_number > 0) {
    # Does any of the insertions affect a gene (where code != 0)
    genes_with_insertion <- sapply(1:insertions_number,
                                   function(i) chromosome[insertions_position[i]])
    genes_with_insertion <- unique(genes_with_insertion[genes_with_insertion != 0])

    # only survive if all genes with insertion are no longer intact
    survival <- all(intactness_vector[genes_with_insertion] < minimum_intactness)

    if (survival == TRUE) {
      # Deleted genes need to get no intactness left
      intactness_vector[genes_with_insertion] <- -9

      for (i in 1:insertions_number) {
        # Check whether there are insertions in intergenic regions and update the `intactness` vector
        intactness_vector <- intergenicMutationEffect(chromosome_coords = chromosome,
                  intactness                  = intactness_vector,
                  cost_of_intergenic_mutation = cost_of_intergenic_mutation,
                  mutation_coord              = insertions_position[i])
      }

      # Insert the new positions into the chromosome
      for (i in 1:insertions_number) {
        ins_position <- insertions_position[i]
        chromosome   <- c(chromosome[1:ins_position],
                          rep(chromosome[ins_position],
                              times = insertions_size[i]),
                              chromosome[(ins_position + 1):length(chromosome)])
      }
    }
  } else {
    stop("Negative number of insertions??")
  }
  return(list(survival          = survival,
              chromosome        = chromosome,
              intactness_vector = intactness_vector))
}


mutateChromosome <- function(chromosome,
                             intactness_vector,
                             cost_of_intergenic_deletion,
                             cost_of_intergenic_insertion,
                             deletion_rate,
                             insertion_rate,
                             minimum_intactness) {
  survival <- FALSE
  stopifnot(!is.null(chromosome))
  # Start with deletions
  # (we divide the process into steps to precent overlapping between insertions and deletions)
  # The following returns a list with:
  # Survival
  # If survival == TRUE, it updates  and returns the size of the chromosome and the
  # intactness vector for the individual
  simulated_deletions <- simulateDeletions(chromosome = chromosome,
    intactness_vector  = intactness_vector,
    cost_of_intergenic_mutation = cost_of_intergenic_deletion,
    mutation_rate      = deletion_rate,
    minimum_intactness = minimum_intactness)

  if (simulated_deletions$survival) {
    # print(paste("Performing insertion to Chromosome", chromosome_n, "at generation", generation))
    stopifnot(!is.null(simulated_deletions$chromosome))
    # Now simulate insertions, but using the updated chromosome and intactness values
    simulated_insertions <- simulateInsertions(chromosome = simulated_deletions$chromosome,
      intactness_vector  = simulated_deletions$intactness_vector,
      cost_of_intergenic_mutation = cost_of_intergenic_insertion,
      mutation_rate      = insertion_rate,
      minimum_intactness = minimum_intactness)

    # If chromosome still survives, add it to the 'germline' of the next generation
    if (simulated_insertions$survival) {
      # Update the record for this chromosome
      survival          <- simulated_insertions$survival
      chromosome        <- simulated_insertions$chromosome
      intactness_vector <- simulated_insertions$intactness_vector
    }
  }
  return(list(survival          = survival,
              chromosome        = chromosome,
              intactness_vector = intactness_vector))
}

# ---------------------------------------------------------------------------- #
# RUN SIMULATION
# ---------------------------------------------------------------------------- #

# Make each chromosome in the population
# Each chromosome is positional vector,
   # where 0 is intergenic region and values > 0 are each gene
# 1 1 0 0 0 0 0 0 0 0 2 2 0 0 0 0 0 0 0 0 3 3 0 0 0 0 0 0 0 0
chromosome        <- lapply(1:gene_number,
  function(i) c(rep(0, intergene_length), rep(i, gene_length)))
chromosome        <- do.call(c, chromosome)

chromosome_list   <- lapply(1:N, function(i) return(chromosome))

# Setup a matrix with the with the 'intactness index'
     # for each gene in each chromosome
# where 1 == 100%, or completely intact
# The genes are numbered 1:gene_number;
# the intactness of the i-th gene is in the i-th position in this vector
intactness_matrix <- t(sapply(1:N, function(i) rep(1, times = gene_number)))


# ---------------------------------------------------------------------------- #
# Save out initial conditions

generation <- 0

print(paste(Sys.time(), ": Simulation for ", generation,
            " finished, writing out.", sep = ''))

chrom_list_out <- paste(chromosome_list_output, generation, sep="/")
chrom_list_out <- paste(chrom_list_out, ".Rda", sep= "")

print(paste(Sys.time(), ": Saving chromosome list to ",
            chrom_list_out, ".", sep = ''))
save(chromosome_list, file = chrom_list_out)

print(paste(Sys.time(), ": Appending intactness matrix to to ",
            intactness_matrix_output, ".", sep = ''))
write.table(intactness_matrix,
            file = intactness_matrix_output,
            append = TRUE,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# ---------------------------------------------------------------------------- #
# Run simulation for each chromosome (i.e. individual) for each generation
for (generation in 1:generation_number) {

  print(paste(Sys.time(), ": Performing generation ", generation, sep = ''))

  # Run simulation for each chromosome in the population
  #  For each chromosome, returns a list with
    # survival (TRUE or FALSE, has the chromosome survived?)
    # new chromosome
    # new intactness vector
   chromosome_result_list <- lapply(1:N,
      function(n) mutateChromosome(chromosome  = chromosome_list[[n]],
                           intactness_vector            = intactness_matrix[n,],
                           cost_of_intergenic_deletion  = cost_of_intergenic_deletion,
                           cost_of_intergenic_insertion = cost_of_intergenic_insertion,
                           deletion_rate                = deletion_rate,
                           insertion_rate               = insertion_rate,
                           minimum_intactness           = minimum_intactness))
  stopifnot(length(chromosome_result_list) == N)
  stopifnot(!sapply(chromosome_result_list, is.null))

  # Which chromosomes survive?
  survivors <- sapply(chromosome_result_list, function(result) result$survival)
  survivors <- c(1:N)[survivors]

  stopifnot(length(survivors) > 0)

  # Select the individuals for the next generation from the survivor pool
  new_generation_chromosomes <- sample(x = survivors, size = N, replace = TRUE)

  chromosome_list   <- lapply(chromosome_result_list[new_generation_chromosomes],
                              function(result) result$chromosome)
  intactness_matrix <- lapply(chromosome_result_list[new_generation_chromosomes],
                              function(result) result$intactness_vector)
  intactness_matrix <- do.call(rbind, intactness_matrix)

  # Write out
  if (generation %in% save_every_nth_seq) {

    print(paste(Sys.time(), ": Simulation for ", generation,
                " finished, writing out.", sep = ''))

    chrom_list_out <- paste(chromosome_list_output, generation, sep="/")
    chrom_list_out <- paste(chrom_list_out, ".Rda", sep= "")

    print(paste(Sys.time(), ": Saving chromosome list to ",
                chrom_list_out, ".", sep = ''))
    save(chromosome_list, file = chrom_list_out)

    print(paste(Sys.time(), ": Appending intactness matrix to to ",
                intactness_matrix_output, ".", sep = ''))
    write.table(intactness_matrix,
                file = intactness_matrix_output,
                append = TRUE,
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE)

  }
}
