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
# Load the functions that run the simulation

source('simulation_functions.r')

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
