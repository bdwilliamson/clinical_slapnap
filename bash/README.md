# The `bash` directory

This directory contains all `bash` code necessary to reproduce the simulation analyses and the SLAPNAP analyses. All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. The simulations were designed to run on a high-performance cluster computer running the Slurm job scheduling system. If you want to run the code locally, or use a different job scheduler, modify the code in the appropriate places to do so (the code is flagged with comments).

# SLAPNAP analyses

To run the SLAPNAP analyses for continuous IC80, run the following:
```{bash}
chmod u+x
./run_slapnap_continuous.sh
```
To run the analyses for binary IC80, run the following:
```{bash}
chmod u+x
./run_slapnap.sh
```
Both bash scripts require that you have Docker installed. There may be additional work to specify the correct Docker arguments if you are running on Windows. See the documentation [here](https://benkeser.github.io/slapnap/) for more information.

# Simulation 1: relative efficiency for estimating the mean value

This simulation requires the following code:
* `run_simulation_1.sh`: submit a batch job array using the given parameters
* `submit_simulation_1.sh`: submit several batch job arrays for each combination of parameters.

Run the simulation using the bash script `submit_simulation_1.sh`, e.g.,
```{bash}
chmod u+x
./submit_simulation_1.sh
```

# Simulation 1b: relative efficiency for estimating the mean value

Run the simulation using the bash script `submit_simulation_1b.sh`, e.g.,
```{bash}
chmod u+x
./submit_simulation_1b.sh
```

# Simulation 2: does SLAPNAP improve sieve analysis?

This simulation requires the following code:
* `run_simulation_2.sh`: submit a batch job array using the given parameters
* `submit_simulation_2.sh`: submit several batch job arrays for each combination of parameters.

Run the simulation using the bash script `submit_simulation_2.sh`, e.g.,
```{bash}
chmod u+x
./submit_simulation_2.sh
```
