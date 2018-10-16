# Simulate chromosome size increase in the absence of recombination

## Running the model

Insertions or deletions in intergenic regions decrease the 'intactness' of adjacent genes, but are not selected against (so they accumulate following drift). Deletions and insertions falling inside genes cause are lethal, but only if the gene is still 30% intact. Deletions are more likely to fall in a gene than an insertions, because their whole length affect the chromosome, whereas insertions only affect the exact position where they fall. Thus, deletions are more likely to be selected against than insertions. If the mutations rate of insertions and deletions is the same, we thus predict that the chromosomes in the population should increase on average until a large number of genes have lost their intactness, where deletions stop being selected against. If deletions and insertions happen at the same rate, we thus predict the size of chromosome to remain constant from this point. If instead deletions happen at a greater frequency than insertions, we expect the chromosome to retract after increasing in size.

We thus run two simulations. In the first, we have an equal mutations rate for insertions and deletions. We select a value that gives us an average of four insertions and four deletions per chromosome per generation.

```sh

mkdir -p tmp
module load R

Rscript simulate.r \ \
	--insertion_rate 5e-6 \ \
	--deletion_rate 5e-6 \ \
	--intactness_matrix tmp/2018-10-14_5e-6-5e-6-matrix.txt \ \
	--chromosome_list_dir tmp/2018-10-14_5e-6-5e-6-matrix-chrom_out \
  1> tmp/2018-10-14_5e-6-5e-6.out \
  2> tmp/2018-10-14_5e-6-5e-6.err

```

In the second simulation, we select a value that gives us double the deletions (eight per chromosome per generation) than insertions (four per chromosome per generation).

```sh

Rscript simulate.r \ \
	--insertion_rate 5e-6 \ \
	--deletion_rate 1e-05 \ \
	--intactness_matrix tmp/2018-10-14_5e-6-1e-05-matrix.txt \ \
	--chromosome_list_dir tmp/2018-10-14_5e-6-1e-05-matrix-chrom_out \
  1> tmp/2018-10-14_5e-6-1e-05.out \
  2> tmp/2018-10-14_5e-6-1e-05.err

```

To increase the run speed, we are only recording the simulation results at every 20th generation.

## Secind simulation

This simulation did not quite work. I have developed a different approach.

## Run simulations

We first test cases with neutral evolution. We care at looking at different rates of insertions and deletions relative to each other.

```sh

mkdir -p tmp/simulations/neutral_1

parallel -j 20 \
'Rscript simulate_2.r \
	--generation_number 500 \
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
	--chromosome_list_out tmp/simulations/neutral_1/simulation_{}.txt' &

	--generation_number 500 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 2/2000 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral TRUE \
	--chromosome_list_out
 \
	--generation_number 40 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1/2000 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 2/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral TRUE \
	--chromosome_list_out

```

We now test whether accruing a cost to deletions increases the size of the chromosome. To do this, we give each locus a fitness value (fixed at 10), which is removed by deletions. Insertions have no fitness cost.

```sh
 \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1/2000 \
	--deletion_size 50 \
	--point_mutation_rate 0 \
	--point_mutation_cost 0 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
	--chromosome_list_out

```

As expected, the average chromosome eventually gets quite large, as the population accumulates insertions but purges highly deleterious deletions.

However, a non-recombining chromosome is thought to undergo sequence-level degeneration as well degeneration based on insertions and deletions. In the following simulation, we include weakly deleterious point mutations that decrease the fitness value of the locus (which starts at 10) by 2. In this simulation, the size of the chromosome eventually plateaus, as the cost of deletions decreases.

```sh
 \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1/2000 \
	--deletion_size 50 \
	--point_mutation_rate 5/200 \
	--point_mutation_cost 2 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
	--chromosome_list_out
 \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1/2000 \
	--deletion_size 50 \
	--point_mutation_rate 5/200 \
	--point_mutation_cost 5 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
	--chromosome_list_out
 \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1/2000 \
	--deletion_size 50 \
	--point_mutation_rate 10/200 \
	--point_mutation_cost 10 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
	--chromosome_list_out

```

Our simulations suggest that the size of the chromosome undergoes a very large increase, but slower to decrease. This is due to the absence of a cost for insertions: after the degeneration of the chromosome, deletions become neutral and the size of the chromosome will evolve mainly via drift. There are two conditions in which the size of the chromosomes can decrease faster after their initial increase. The first a having a deletion rate that is slightly higher than the insertion rate, which slows down the initial increase in size of the non-recombining chromosome, but eventually accelerates its decrease. The other is having cost associated with 'junk' loci in the genome. We simulated this by adding a multiplier cost that reduces an individual's fitness proportionally to the number of non-functional loci (i.e. the loci introduced by insertions or rendered non-functional by repeated point mutations).


```sh

 \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1/2000 \
	--deletion_size 50 \
	--point_mutation_rate 5/200 \
	--point_mutation_cost 5 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
	--chromosome_list_out
 \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1.2/2000 \
	--deletion_size 50 \
	--point_mutation_rate 5/200 \
	--point_mutation_cost 5 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0 \
	--every_nth 5 \
	--neutral FALSE \
	--chromosome_list_out
 \
	--generation_number 1000 \
	--population_size 1000 \
	--chromosome_length 4000 \
	--deletion_rate 1/2000 \
	--deletion_size 50 \
	--point_mutation_rate 5/200 \
	--point_mutation_cost 5 \
	--locus_value 10 \
	--insertion_rate 1/2000 \
	--insertion_size 50 \
	--insertion_cost 0.5 \
	--every_nth 5 \
	--neutral FALSE \
	--chromosome_list_out

```
