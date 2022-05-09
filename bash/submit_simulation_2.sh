#!/bin/bash

# run the entire simulation 1
ml fhR/4.0.2-foss-2019b
ml jbigkit

save_dir="/fh/fast/gilbert_p/user/bwillia2/clinical_slapnap/simulation_2/"
mkdir -p $save_dir
nreps_per_job=100

analyses=("priority" "site-scanning")
for analysis in ${analyses[@]}; do
  io_prefix="clinical_slapnap/output_${analysis}"
  ./run_simulation_2.sh $analysis 1000 $nreps_per_job $save_dir $io_prefix
done
