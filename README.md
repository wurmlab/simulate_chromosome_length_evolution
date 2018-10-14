# Simulate chromosome size increase in the absence of recombination

## Running the model

Insertions or deletions in intergenic regions decrease the 'intactness' of adjacent genes, but are not selected against (so they accumulate following drift). Deletions and insertions falling inside genes cause are lethal, but only if the gene is still 30% intact. Deletions are more likely to fall in a gene than an insertions, because their whole length affect the chromosome, whereas insertions only affect the exact position where they fall. Thus, deletions are more likely to be selected against than insertions. If the mutations rate of insertions and deletions is the same, we thus predict that the chromosomes in the population should increase on average until a large number of genes have lost their intactness, where deletions stop being selected against. If deletions and insertions happen at the same rate, we thus predict the size of chromosome to remain constant from this point. If instead deletions happen at a greater frequency than insertions, we expect the chromosome to retract after increasing in size.

We thus run two simulations. In the first, we have an equal mutations rate for insertions and deletions. We select a value that gives us an average of four insertions and four deletions per chromosome per generation.

```sh

mkdir -p tmp
module load R

Rscript simulate.r \
  --insertion_rate 5e-6 \
  --deletion_rate 5e-6 \
  --intactness_matrix tmp/2018-10-14_5e-6-5e-6-matrix.txt \
  --chromosome_list_dir tmp/2018-10-14_5e-6-5e-6-matrix-chrom_out \
  1> tmp/2018-10-14_5e-6-5e-6.out \
  2> tmp/2018-10-14_5e-6-5e-6.err

```

In the second simulation, we select a value that gives us double the deletions (eight per chromosome per generation) than insertions (four per chromosome per generation).

```sh

Rscript simulate.r \
  --insertion_rate 5e-6 \
  --deletion_rate 1e-05 \
  --intactness_matrix tmp/2018-10-14_5e-6-1e-05-matrix.txt \
  --chromosome_list_dir tmp/2018-10-14_5e-6-1e-05-matrix-chrom_out \
  1> tmp/2018-10-14_5e-6-1e-05.out \
  2> tmp/2018-10-14_5e-6-1e-05.err

```

To increase the run speed, we are only recording the simulation results at every 20th generation.
