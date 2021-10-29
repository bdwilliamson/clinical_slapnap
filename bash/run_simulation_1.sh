#!/bin/bash

# run simulation 1

# takes X command-line arguments
# 01: bnab of interest, e.g., "VRC01"
# 02: outcome of interest, e.g., "ic80"
# 03: number of total replicates, e.g., 1000
# 04: number of replicates per job, e.g., 1000
# 05: output directory, e.g., "/fh/fast/gilbert_p/bwillia2/clinical_slapnap/simulation_1/"
# 06: the i/o file
io_prefix="/fh/scratch/delete90/gilbert_p/bwillia2/${6}"
mkdir -p $io_prefix
io_file="$io_prefix/slurm-%A_%a.out"

# run the R script
echo -e \
    '#!/bin/bash\n Rscript simulation_1.R ' \
    '--bnab $1 --outcome $2 --nreps-total $3 ' \
    '--nreps-per-job $4 --output-dir $5' > run_sim1.sh
chmod u+x run_sim1.sh

njobs=`expr $3 / $4 \* 3`
arry="1-$njobs"

sbatch -A gilbert_p --time=1-0 --array=$arry -e $io_file \
    -o $io_file ./run_sim1.sh $1 $2 $3 $4 $5
rm run_sim1.sh
