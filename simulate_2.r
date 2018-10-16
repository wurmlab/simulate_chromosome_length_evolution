#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------- #
# DESCRIPTION HERE



# ---------------------------------------------------------------------------- #

library("optparse")

option_list <- list(
  make_option(c("-g", "--generation_number"), type="numeric", default=20,
		          help="Number of generations [default %default]", metavar="number"),
  make_option(c("-N", "--population_size"), type="numeric", default=100,
		          help="Number of haploid individuals (chromosomes) [default %default]", metavar="number"),
  make_option(c("-L", "--chromosome_length"), type="numeric", default=1000,
		          help="Chromosome length [default %default]", metavar="number"),
  make_option(c("-l", "--locus_value"), type="numeric", default=10,
		          help="Value of each locus in the chromosome [default %default]", metavar="number"),
  make_option(c("-d", "--deletion_rate"), type="numeric", default=0.001,
		          help="Deletion rate [default %default]", metavar="number"),
  make_option(c("-D", "--deletion_size"), type="numeric", default=10,
		          help="Deletion size [default %default]", metavar="number"),
  make_option(c("-i", "--insertion_rate"), type="numeric", default=0.001,
		          help="Insertion rate [default %default]", metavar="number"),
  make_option(c("-I", "--insertion_size"), type="numeric", default=10,
		          help="Insertion size [default %default]", metavar="number"),
  make_option(c("-j", "--insertion_cost"), type="numeric", default=10,
		          help="Cost of each base-pair of insertion [default %default]", metavar="number"),
  make_option(c("-p", "--point_mutation_rate"), type="numeric", default=1,
		          help="Insertion size [default %default]", metavar="number"),
  make_option(c("-k", "--point_mutation_cost"), type="numeric", default=10,
		          help="Cost of each point mutations [default %default]", metavar="number"),
  make_option("--neutral", type="logical", default=FALSE,
		          help="Remove the costs of all mutations [default %default]", metavar="number"),
  make_option("--every_nth", type="numeric", default=1,
		          help="Save only every n-th simulation [default %default]", metavar="number"),
  make_option(c("-o", "--chromosome_list_out"), type="character", default="out.Rda",
		          help="File for output chromosome list [default= %default]", metavar="character")
)

# ---------------------------------------------------------------------------- #

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #

simulateDeletions <- function(chromosome,
                              mutation_rate,
                              mutation_size,
                              chromosome_length) {

  if (chromosome_length > mutation_size) {
    # Set up the mutations
    ## generate the number of mutations
    mutation_number   <- rpois(n = 1,
                               lambda = mutation_rate * chromosome_length)

    if (mutation_number > 0) {
      ## Mutation positions
      deletion_starts <- sample(x       = 1:chromosome_length,
                                size    = mutation_number,
                                replace = FALSE)

      deletion_ends   <- deletion_starts + mutation_size

      for (i in 1:mutation_number) {
        chromosome <- chromosome[-c(deletion_starts[i]:deletion_ends[i])]
      }
    }
  }
  return(chromosome)
}

simulatePointMutations <- function(chromosome,
                                   mutation_rate,
                                   mutation_cost,
                                   original_chromosome_length) {
  # Set up the mutations
  if (mutation_rate > 0) {
    ## generate the number of mutations
    mutation_number <- round(mutation_rate * original_chromosome_length)

    if (mutation_number > 0) {
      ## Mutation positions
      mutation_position <- sample(x       = 1:length(chromosome),
                                  size    = mutation_number,
                                  replace = FALSE)
      # Reduce value of mutated sites to 0, unless it has already been reduced by
          # other means
      site_values                   <- chromosome[mutation_position]
      site_values[site_values > 0]  <- site_values[site_values > 0] - mutation_cost
      site_values[site_values < 0]  <- 0
      chromosome[mutation_position] <- site_values
    }
  }
  return(chromosome)
}

simulateInsertions <- function(chromosome,
                               mutation_rate,
                               mutation_size,
                               original_chromosome_length) {

  # Set up the mutations
  ## generate the number of mutations
  # Edge effect: cannot have an insertion in the last locus
  mutation_number <- rpois(n = 1, lambda = mutation_rate * original_chromosome_length)

  if (mutation_number > 0) {

    ## Mutation positions
    mutation_starts <- sample(x       = 1:(length(chromosome) - 1),
                              size    = mutation_number,
                              replace = FALSE)
    # For each insertion, slice the chromosome at position, add insertions,
          # and stich back together
    for (i in 1:mutation_number) {
      ins_position <- mutation_starts[i]
      chromosome   <- c(chromosome[1:ins_position],
                        rep(0, times = mutation_size),
                        chromosome[(ins_position + 1):length(chromosome)])
    }
  }
  return(chromosome)
}

simulateMutations <- function(chromosome,
                              deletion_rate,
                              deletion_size,
                              point_mutation_rate,
                              point_mutation_cost,
                              insertion_rate,
                              insertion_size) {

  chromosome_length <- length(chromosome)
  if (chromosome_length > 0) {
    # Deletions, then point mutations, then insertions
    if (deletion_rate > 0) {
          chromosome <- simulateDeletions(chromosome        = chromosome,
                                          mutation_rate     = deletion_rate,
                                          mutation_size     = deletion_size,
                                          chromosome_length = chromosome_length)
    }

    if (point_mutation_rate > 0) {
      chromosome <- simulatePointMutations(chromosome        = chromosome,
                                           mutation_rate     = point_mutation_rate,
                                           mutation_cost     = point_mutation_cost,
                                           original_chromosome_length = chromosome_length)
    }

    if (insertion_rate > 0) {
      chromosome <- simulateInsertions(chromosome        = chromosome,
                                       mutation_rate     = insertion_rate,
                                       mutation_size     = insertion_size,
                                       original_chromosome_length = chromosome_length)
    }
  }
  return(chromosome)
}

runSimulation <- function(generation_number,
                          N,
                          chromosome_length,
                          deletion_rate,
                          deletion_size,
                          point_mutation_rate,
                          point_mutation_cost,
                          locus_value,
                          insertion_rate,
                          insertion_size,
                          insertion_cost,
                          every_nth,
                          neutral_mode = FALSE) {

  # Make each chromosome in the population
  # Each chromosome is a vector of genetic values, with a base level of 1
  # 1 1 1 1 1 1 1 1 1 1 1 1 1 1
  # Substitutions decrease the number to 0
  # Deletions remove the genetic values
  # Insertions add values (0 if neutral, negative number if costly)
  chromosome        <- rep(locus_value, times = chromosome_length)
  chromosome_list   <- lapply(1:N, function(i) return(chromosome))

  results_list      <- list()
  results_list[[1]] <- data.frame(generation = 1,
                                  chromLenSummary(chromosome_list, locus_value))

  save_every_nth_seq <- seq(from = every_nth, to = generation_number, by = every_nth)

  for (generation in 1:generation_number) {

    stopifnot(insertion_cost >= 0 & insertion_cost <= 1)

    # print(paste(Sys.time(), ": Performing generation ", generation, sep = ''))

    # Step one: mutation
    simulation_results    <- lapply(chromosome_list, function(chromosome)
       simulateMutations(chromosome          = chromosome,
                         deletion_rate       = deletion_rate,
                         deletion_size       = deletion_size,
                         point_mutation_rate = point_mutation_rate,
                         point_mutation_cost = point_mutation_cost,
                         insertion_rate      = insertion_rate,
                         insertion_size      = insertion_size))

    stopifnot(!sapply(simulation_results, function(x) any(is.na(x))))

    if (any(sapply(simulation_results, length) > (5 * chromosome_length))) {
      warning("Some chromosomes too large, breaking early")
      break
    }

    # Step two: selection
      # Sum the fitness of each locus for each individual
      # Scale it relative to the individual with highest fitness
      # Select N individuals with probability based on those fitnesses
    if (neutral_mode == FALSE) {
      chromosome_fitness <- sapply(simulation_results, sum)

      # Cost of insertions
       # It's a multiplier cost - if 40% of the chromosome is '0',
        # and the cost is 1, then we reduce the fitness of the individual by 40%
      if (insertion_cost > 0) {
        zero_length                 <- sapply(simulation_results,
                                              function(x) sum(x == 0)/length(x))
        insertion_fitness_reduction <- zero_length * insertion_cost
        chromosome_fitness          <- chromosome_fitness - insertion_fitness_reduction * chromosome_fitness
      }

      chromosome_fitness <- chromosome_fitness/max(chromosome_fitness)

    } else if (neutral_mode == TRUE) {
      chromosome_fitness <- rep(1/length(simulation_results), times = length(simulation_results))
    } else {
      stop("neutral_mode can only be TRUE or FALSE")
    }

    stopifnot(chromosome_fitness >= 0 & chromosome_fitness <= 1)


    next_generation <- sample(x       = 1:N,
                              size    = N,
                              prob    = chromosome_fitness,
                              replace = TRUE)

    chromosome_list <- simulation_results[next_generation]

    if (generation %in% save_every_nth_seq) {
      generation_result_df            <- chromLenSummary(chromosome_list, locus_value)
      generation_result_df$generation <- generation
      col_order <- colnames(generation_result_df)
      col_order <- c('generation', col_order[col_order != 'generation'])
      generation_result_df <- generation_result_df[, col_order]
      results_list[[which(save_every_nth_seq == generation) + 1]] <- generation_result_df
    }
  }

  results_list <- do.call(rbind, results_list)
  return(results_list)

}

chromLenSummary <- function(chromosome_list, locus_value) {

  chr_lengths      <- sapply(chromosome_list, length)

  chr_value <- sapply(chromosome_list, sum)
  chr_junk  <- sapply(chromosome_list, function(x) sum(x == 0))
  chr_delet <- sapply(chromosome_list, function(x) sum(x < locus_value))


  return(data.frame(
    mean_length   = mean(chr_lengths),
    median_length = median(chr_lengths),
    sd_length     = sd(chr_lengths),
    mad_length    = mad(chr_lengths),
    max_length    = max(chr_lengths),
    min_length    = min(chr_lengths),
    mean_value    = mean(chr_value),
    sd_vale       = sd(chr_value),
    median_value  = median(chr_value),
    mad_value     = mad(chr_value),
    max_value     = max(chr_value),
    min_value     = min(chr_value),
    mean_junk     = mean(chr_junk),
    sd_junk       = sd(chr_junk),
    median_junk   = median(chr_junk),
    mean_mutated_loci = mean(chr_delet),
    median_mutated_loci = median(chr_delet),
    sd_mutated_loci     = median(chr_delet)
  ))
}

# ---------------------------------------------------------------------------- #
# RUN SIMULATION
# ---------------------------------------------------------------------------- #

simulation <- runSimulation(generation_number   = opt$generation_number,
                                 N                   = opt$population_size,
                                 chromosome_length   = opt$chromosome_length,
                                 deletion_rate       = opt$deletion_rate,
                                 deletion_size       = opt$deletion_size,
                                 point_mutation_rate = opt$point_mutation_rate,
                                 point_mutation_cost = opt$point_mutation_cost,
                                 locus_value         = opt$locus_value,
                                 insertion_rate      = opt$insertion_rate,
                                 insertion_size      = opt$insertion_size,
                                 insertion_cost      = opt$insertion_cost,
                                 every_nth           = opt$every_nth,
                                 neutral             = opt$neutral)

# ---------------------------------------------------------------------------- #

write.table(simulation,
            file      = opt$chromosome_list_out,
            quote     = FALSE,
            col.names = FALSE,
            row.names = FALSE)

# ---------------------------------------------------------------------------- #
