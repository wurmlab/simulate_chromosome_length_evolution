# Simulate chromosome size evolution in the absence of recombination

Using optical maps and whole-genome sequences, we have shown that the supergene variant with suppressed recombination has increased in size. Our interpretations of these results is that this supergene variant has undergone "degenerative expansion". Our argument is that although non-recombining chromosomes accumulate deleterious mutations, they are expected to accumulate weak deleterious mutations faster than relatively stronger deleterious mutations. Because deletions are thought to be generally more deleterious than insertions, non-recombining chromosomes are expected to increase in size.

We have produced a set of simulations to illustrate our argument. Briefly, the simulations include a population of 1000 individuals, each carrying a non-recombining chromosome composed of "functional" elements. These chromosomes are affected by insertions and deletions of the same size, at the same rate. Deletions are deleterious because they remove functional elements. We can give an arbitrary cost to insertions, as well as background degeneration and an additional class of very large deletions. The simulations show that the the average chromosome size in a population increases where deletions are more deleterious than insertions.

## Run simulations

The functions used for the simulations are in the script `simulation_functions.r`. Each simulation is called with the script `run_simulation.r`, which includes the option of modifying many parameters. The script gets called in `simulation_list.sh` for each simulation set, which uses GNU-parallel. This script also groups similar simulation sets into single files, to then load in R for further analysis.

```sh

mkdir -p tmp/simulations
mkdir -p results/fig

sh simulation_list.sh \
	1> tmp/simulation_list.sh.out \
	tmp/simulation_list.sh.err

```

## Analysis of simulations

The analysis of the simulations is done in knitr in the file `simulation_report.Rmd`, which produced the file `simulation_report.html`.
