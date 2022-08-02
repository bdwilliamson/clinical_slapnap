#!/bin/bash

# run simulation 2

# takes X command-line arguments
# 01: the analysis (priority, site-scanning)
# 02: number of total replicates, e.g., 1000
# 03: number of replicates per job, e.g., 1000
# 04: output directory, e.g., "/fh/fast/gilbert_p/bwillia2/clinical_slapnap/simulation_1/"
# 05: i/o prefix
io_prefix="/fh/scratch/delete90/gilbert_p/bwillia2/${5}"
mkdir -p $io_prefix
io_file="$io_prefix/slurm-%A_%a.out"

# run the R script
echo -e \
    '#!/bin/bash\n Rscript simulation_2.R ' \
    '--analysis $1 --nreps-total $2 ' \
    '--nreps-per-job $3 --output-dir $4' > run_sim2.sh
chmod u+x run_sim2.sh

njobs=`expr $2 / $3 \* 12`
arry="1-$njobs"

# modify the following line to run on your cluster environment
sbatch -A gilbert_p --time=1-0 --array=$arry -e $io_file \
    -o $io_file ./run_sim2.sh $1 $2 $3 $4
rm run_sim2.sh
