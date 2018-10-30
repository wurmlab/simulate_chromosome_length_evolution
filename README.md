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

## Second simulation

This simulation did not quite work. I have developed a different approach.

## Run simulations

We first test cases with neutral evolution. We care at looking at different rates of insertions and deletions relative to each other.

```sh

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

```

We now test whether accruing a cost to deletions increases the size of the chromosome. To do this, we give each locus a fitness value (fixed at 10), which is removed by deletions. Insertions have no fitness cost.

```sh

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


```

As expected, the average chromosome eventually gets quite large, as the population accumulates insertions but purges highly deleterious deletions.

## Background degeneration

However, a non-recombining chromosome is thought to undergo sequence-level degeneration as well degeneration based on insertions and deletions. In the following simulation, we include weakly deleterious point mutations that decrease the fitness value of the locus (which starts at 10) by 2. In this simulation, the size of the chromosome eventually plateaus, as the cost of deletions decreases.

```sh

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

```

Our simulations suggest that the size of the chromosome undergoes a very large increase, but slower to decrease. This is due to the absence of a cost for insertions: after the degeneration of the chromosome, deletions become neutral and the size of the chromosome will evolve mainly via drift. There are two conditions in which the size of the chromosomes can decrease faster after their initial increase. The first a having a deletion rate that is slightly higher than the insertion rate, which slows down the initial increase in size of the non-recombining chromosome, but eventually accelerates its decrease. The other is having cost associated with 'junk' loci in the genome. We simulated this by adding a multiplier cost that reduces an individual's fitness proportionally to the number of non-functional loci (i.e. the loci introduced by insertions or rendered non-functional by repeated point mutations).


```sh

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
	awk '{print "deletion_bias_no_background_degeneration2\t"FILENAME"\t"$0}' $p >> ../insertion_cost_all_sim
done
cd ../../../

cd tmp/simulations/deletion_bias_no_background_degeneration1
for p in *; do
	awk '{print "deletion_bias_no_background_degeneration1\t"FILENAME"\t"$0}' $p >> ../insertion_cost_all_sim
done
cd ../../../

```

## Effect of having a few very large deletions

```sh

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

```
