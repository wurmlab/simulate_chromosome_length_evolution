
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
                               insertion_cost,
                               original_chromosome_length) {
  stopifnot(length(chromosome) > 4)
  # Set up the mutations
  ## generate the number of mutations
  # Edge effect: cannot have an insertion in the last locus
  mutation_number <- rpois(n = 1, lambda = mutation_rate * original_chromosome_length)

  ## edge case: more inserts than base pairs
  if (mutation_number > (length(chromosome) - 1)) {
    mutation_number = mutation_number - 1
  }

  if (mutation_number > 0) {

    ## Mutation positions
        ## Because of edge effects, we cannot have insertions in the first or
        ## last two loci in the chromosome
    mutation_starts <- sample(x       = 2:(length(chromosome) - 2),
                              size    = mutation_number,
                              replace = FALSE)
    mutation_starts <- sort(mutation_starts)
    # For each insertion, slice the chromosome at position, add insertions,
          # and stich back together
    for (i in 1:mutation_number) {
      ins_position <- mutation_starts[i]
      # insertion cost affects the locus just before and just after the insertions
      start_value <- chromosome[ins_position] - insertion_cost
      start_value <- ifelse(start_value < 0, 0, start_value)

      end_value   <- chromosome[ins_position + 1] - insertion_cost
      end_value   <- ifelse(end_value < 0, 0, end_value)

      # Add the loci with cost and the insertion
      chromosome   <- c(chromosome[1:(ins_position - 1)],
                        start_value,
                        rep(0, times = mutation_size),
                        end_value,
                        chromosome[(ins_position + 2):length(chromosome)])

      # update mutation_starts with the length of the insertion
      if (i < mutation_number) {
        i_remaining <- (i+1):length(mutation_starts)
        mutation_starts[i_remaining] <- mutation_starts[i_remaining] + mutation_size
      }
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
                              insertion_size,
                              insertion_cost) {

  chromosome_length <- length(chromosome)
  # Deletions, then point mutations, then insertions
  # We also need to be careful with edge cases --- inserting does not work in
      # chromosomes with  single base pair ---
      # if a chromosome has a single base pair, force it to NULL
  if (chromosome_length > 0) {
    if (deletion_rate > 0) {
          chromosome <- simulateDeletions(chromosome        = chromosome,
                                          mutation_rate     = deletion_rate,
                                          mutation_size     = deletion_size,
                                          chromosome_length = chromosome_length)
    }
  }

  # Insertions do not work in chromosomes with smaller than 5 loci ---
      # if a chromosome has a 4 loci or less, force it to 0 loci
  if (length(chromosome) <= 4) {
    chromosome <- chromosome[-1]
  }

  if (length(chromosome) > 0) {
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
                                       insertion_cost    = insertion_cost,
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
                          junk_cost,
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

    stopifnot(junk_cost >= 0 & junk_cost <= 1)

    # print(paste(Sys.time(), ": Performing generation ", generation, sep = ''))

    # Step one: mutation
    simulation_results    <- lapply(chromosome_list, function(chromosome)
       simulateMutations(chromosome          = chromosome,
                         deletion_rate       = deletion_rate,
                         deletion_size       = deletion_size,
                         point_mutation_rate = point_mutation_rate,
                         point_mutation_cost = point_mutation_cost,
                         insertion_rate      = insertion_rate,
                         insertion_size      = insertion_size,
                         insertion_cost      = insertion_cost))

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
      if (junk_cost > 0) {
        zero_length                 <- sapply(simulation_results,
                  function(x) ifelse(length(x) > 0, sum(x == 0)/length(x), 0))
        insertion_fitness_reduction <- zero_length * junk_cost
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
