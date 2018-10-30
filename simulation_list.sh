#!/usr/bin/env bash

# We first test cases with neutral evolution.
# We care at looking at different rates of insertions and deletions relative to each other.

module load parallel/20170422
module load R/3.4.3

mkdir -p tmp/simulations/neutral_1
parallel -j 25 \
'Rscript run_simulations.r \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral TRUE \
	--chromosome_list_out tmp/simulations/neutral_1/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/neutral_2
parallel -j 25 \
'Rscript run_simulations.r \
	--generation_number 500 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 0.001 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral TRUE \
  --chromosome_list_out tmp/simulations/neutral_2/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/neutral_3
parallel -j 25 \
'Rscript run_simulations.r \
	--generation_number 40 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 0.001 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral TRUE \
  --chromosome_list_out tmp/simulations/neutral_3/simulation_{}.txt' ::: {1..25}

## HIgher cost for deletions than insertions

# We now test whether accruing a cost to deletions increases the size of the chromosome.
# To do this, we give each locus a fitness value (fixed at 10), which is removed by deletions.
# Insertions have no fitness cost.

 mkdir -p tmp/simulations/deletion_with_cost
 parallel -j 25 \
 'Rscript run_simulations.r \
 	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_with_cost/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/deletion_with_cost_insertions_with_cost
parallel -j 25 \
'Rscript run_simulations.r \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--insertion_cost 10 \
	--every_nth 5 \
	--neutral FALSE \
 --chromosome_list_out tmp/simulations/deletion_with_cost_insertions_with_cost/simulation_{}.txt' ::: {1..25}

## Backgroung degeneration

# a non-recombining chromosome is thought to undergo sequence-level degeneration
# as well degeneration based on insertions and deletions. In the following
# simulation, we include weakly deleterious point mutations that decrease the
# fitness value of the locus (which starts at 10) by 2. In this simulation, the
# size of the chromosome eventually plateaus, as the cost of deletions decreases.

mkdir -p tmp/simulations/point_mutations_2
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0.025 \
	--point_mutation_cost 5 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
  --chromosome_list_out tmp/simulations/point_mutations_2/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/point_mutations_3
parallel -j 25 \
'Rscript run_simulations.r \
  --generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0.05 \
	--point_mutation_cost 10 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
  --chromosome_list_out tmp/simulations/point_mutations_3/simulation_{}.txt' ::: {1..25}

# Our simulations suggest that the size of the chromosome undergoes a very large
#  increase, but slower to decrease. This is due to the absence of a cost for
# insertions: after the degeneration of the chromosome, deletions become neutral
# and the size of the chromosome will evolve mainly via drift.
# There are two conditions in which the size of the chromosomes can decrease
# faster after their initial increase. The first a having a deletion rate that
# is slightly higher than the insertion rate, which slows down the initial
# increase in size of the non-recombining chromosome, but eventually accelerates
# its decrease. The other is having cost associated with 'junk' loci in the
# genome. We simulated this by adding a multiplier cost that reduces an
# individual's fitness proportionally to the number of non-functional loci
# (i.e. the loci introduced by insertions or rendered non-functional by
# repeated point mutations).

mkdir -p tmp/simulations/deletion_bias
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 6e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0.025 \
  --point_mutation_cost 5 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_bias/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/deletion_bias_no_background_degeneration
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 6e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_bias_no_background_degeneration/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/deletion_bias_no_background_degeneration1
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5.5e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_bias_no_background_degeneration1/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/deletion_bias_no_background_degeneration2
parallel -j 15 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5.1e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_bias_no_background_degeneration2/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/deletion_size_bias_no_background_degeneration
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5e-04 \
  --deletion_size 75 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_size_bias_no_background_degeneration/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/deletion_bias_insertion_cost_no_background_degeneration
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 6e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 10 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_bias_insertion_cost_no_background_degeneration/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/deletion_bias_insertion_cost_background_degeneration
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 6e-04 \
  --deletion_size 50 \
	--point_mutation_rate 0.025 \
  --point_mutation_cost 5 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 10 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/deletion_bias_insertion_cost_background_degeneration/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/large_deletion_1
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --large_deletion_rate 1e-06 \
  --large_deletion_size 2000 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/large_deletion_1/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/large_deletion_2
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1500 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --large_deletion_rate 2.5e-06 \
  --large_deletion_size 2000 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/large_deletion_2/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/large_deletion_3
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --large_deletion_rate 2.5e-06 \
  --large_deletion_size 5000 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/large_deletion_3/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/large_deletion_4
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1500 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --large_deletion_rate 1.25e-5 \
  --large_deletion_size 2000 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/large_deletion_4/simulation_{}.txt' ::: {1..25}


mkdir -p tmp/simulations/large_deletion_5
parallel -j 25 \
  'Rscript run_simulations.r \
  --generation_number 1000 \
  --population_size 1000 \
  --chromosome_length 4000 \
  --deletion_rate 5e-04 \
  --deletion_size 50 \
  --point_mutation_rate 0 \
  --point_mutation_cost 0 \
  --locus_value 10 \
  --insertion_rate 5e-04 \
  --insertion_size 50 \
  --insertion_cost 0 \
  --large_deletion_rate 1e-06 \
  --large_deletion_size 10000 \
  --every_nth 5 \
  --neutral FALSE \
  --chromosome_list_out tmp/simulations/large_deletion_5/simulation_{}.txt' ::: {1..25}

mkdir -p tmp/simulations/insertion_cost_1
parallel -j 25 \
  'Rscript run_simulations.r \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0.025 \
	--point_mutation_cost 5 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--junk_cost 0.5 \
	--every_nth 5 \
	--neutral FALSE \
  --chromosome_list_out tmp/simulations/insertion_cost_1/simulation_{}.txt' ::: {1..25}


mkdir -p tmp/simulations/insertion_cost_2
parallel -j 25 \
  'Rscript run_simulations.r \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 5e-04 \
	--deletion_size 50 \
	--point_mutation_rate 0.025 \
	--point_mutation_cost 5 \
	--locus_value 10 \
	--insertion_rate 5e-04 \
	--insertion_size 50 \
	--junk_cost 1 \
	--every_nth 5 \
	--neutral FALSE \
  --chromosome_list_out tmp/simulations/insertion_cost_2/simulation_{}.txt' ::: {1..25}


# ---------------------------------------------------------------------------- #

# Join up similar simulation sets into single files,
#   to help with following step of loading in R

cd tmp/simulations/neutral_1
for p in *; do
	awk '{print "neutral_1\t"FILENAME"\t"$0}' $p >> ../neutral_all_sim
done
cd ../../../

cd tmp/simulations/neutral_2
for p in *; do
	awk '{print "neutral_2\t"FILENAME"\t"$0}' $p >> ../neutral_all_sim
done
cd ../../../

cd tmp/simulations/neutral_3
for p in *; do
	awk '{print "neutral_3\t"FILENAME"\t"$0}' $p >> ../neutral_all_sim
done
cd ../../../

cd tmp/simulations/deletion_with_cost
for p in *; do
	awk '{print "deletion_with_cost\t"FILENAME"\t"$0}' $p >> ../deletion_with_cost_all_sim
done
cd ../../../

cd tmp/simulations/deletion_with_cost_insertions_with_cost
for p in *; do
	awk '{print "deletion_cost_insertion_cost\t"FILENAME"\t"$0}' $p >> ../deletion_with_cost_all_sim
done
cd ../../../

## Background degeneration

cd tmp/simulations/point_mutations_2
for p in *; do
	awk '{print "point_mutations_2\t"FILENAME"\t"$0}' $p >> ../point_mutations_all_sim
done
cd ../../../

cd tmp/simulations/point_mutations_3
for p in *; do
	awk '{print "point_mutations_3\t"FILENAME"\t"$0}' $p >> ../point_mutations_all_sim
done
cd ../../../

cd tmp/simulations/deletion_bias
for p in *; do
	awk '{print "deletion_bias\t"FILENAME"\t"$0}' $p >> ../insertion_cost_all_sim
done
cd ../../../

cd tmp/simulations/deletion_bias_no_background_degeneration
for p in *; do
	awk '{print "deletion_bias_no_background\t"FILENAME"\t"$0}' $p >> ../insertion_cost_all_sim
done
cd ../../../

cd tmp/simulations/insertion_cost_1
for p in *; do
	awk '{print "insertion_cost_1\t"FILENAME"\t"$0}' $p >> ../insertion_cost_all_sim
done
cd ../../../

cd tmp/simulations/insertion_cost_2
for p in *; do
	awk '{print "insertion_cost_2\t"FILENAME"\t"$0}' $p >> ../insertion_cost_all_sim
done
cd ../../../

cd tmp/simulations/deletion_bias_no_background_degeneration2
for p in *; do
	awk '{print "deletion_bias_no_background_degeneration2\t"FILENAME"\t"$0}' $p \
	  >> ../insertion_cost_all_sim
done
cd ../../../

cd tmp/simulations/deletion_bias_no_background_degeneration1
for p in *; do
	awk '{print "deletion_bias_no_background_degeneration1\t"FILENAME"\t"$0}' $p \
	  >> ../insertion_cost_all_sim
done
cd ../../../


## Effect of having a few very large deletions

cd tmp/simulations/large_deletion_1
for p in *; do
	awk '{print "large_deletion_1\t"FILENAME"\t"$0}' $p >> ../large_deletion_all_sim
done
cd ../../../

cd tmp/simulations/large_deletion_2
for p in *; do
	awk '{print "large_deletion_2\t"FILENAME"\t"$0}' $p >> ../large_deletion_all_sim
done
cd ../../../

cd tmp/simulations/large_deletion_3
for p in *; do
	awk '{print "large_deletion_3\t"FILENAME"\t"$0}' $p >> ../large_deletion_all_sim
done
cd ../../../

cd tmp/simulations/large_deletion_4
for p in *; do
	awk '{print "large_deletion_4\t"FILENAME"\t"$0}' $p >> ../large_deletion_all_sim
done
cd ../../../

cd tmp/simulations/large_deletion_5
for p in *; do
	awk '{print "large_deletion_5\t"FILENAME"\t"$0}' $p >> ../large_deletion_all_sim
done
cd ../../../
