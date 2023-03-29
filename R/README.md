# The `R` directory

This directory contains all `R` code necessary to reproduce the simulation analyses. All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. The simulations were designed to run on a high-performance cluster computer running the Slurm job scheduling system. If you want to run the code locally, or use a different job scheduler, modify the code in the appropriate places to do so (the code is flagged with comments).

# Table 1: summary statistics

Run
```{bash}
Rscript table_1_statistics.R
```
to obtain the summary statistics in Table 1. The code relies on several utility functions that can be found in `00_utils.R`, `02_multi_ab.Rlib`, and `compute_table_1_stats_nab.R`.

# Simulation 1: relative efficiency for estimating the mean value

To reproduce this simulation, the following files are required:
* `summary_stats_for_sim1.R`: create a table of summary statistics necessary for simulation 1.
* `00_sim_utils.R`: functions for generating data and estimating the mean value with auxiliary information.
* `simulation_1.R`: the main script, runs the simulation a specified number of times.

First, run `summary_stats_for_sim1.R`, using, e.g.,
```{bash}
Rscript summary_stats_for_sim1.R
```

Next, run the simulation using the bash script `submit_simulation_1.sh`. Once the analysis is finished, run `load_sim_1.R` to create the plots.

# Simulation 1b: relative efficiency for estimating the mean value

To reproduce this analysis, the following files are required:
* `00_utils.R`: functions useful for SLAPNAP.
* `00_sim_utils.R`: functions for generating data and estimating the mean value with auxiliary information.
* `03_super_learner_libraries.R`: functions for predicting IC80 values.
* `simulation_1b.R`: the main script, runs the simulation a specified number of times.

This analysis requires that you have first run SLAPNAP (specified above) for each bnAb of interest, and the results are saved in `docker_output`.

Then, the analysis is run using the bash script `submit_simulation_1b.sh`. Once the analysis is finished, run `load_sim_1b.R` to create the plots.


# Simulation 2: does SLAPNAP improve sieve analysis?

To reproduce this analysis, the following files are required:
* `00_sim_utils.R`: functions for generating data and estimating the mean value with auxiliary information.
* `lunnMcneil.R`: functions for the Lunn & McNeil test (provided by Michal Juraska).

The analysis is run using the bash script `submit_simulation_2.sh`. Once the analysis is finished, run `load_sim_2.R` to create the plots.
