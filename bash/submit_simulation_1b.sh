#!/bin/bash

# run the entire simulation 1b
ml fhR/4.1.2-foss-2020b
ml jbigkit

# all bnabs
bnabs=("VRC01" "VRC07-523-LS" "PGT121" "VRC26.25" "PGDM1400" \
        "VRC07-523-LS_PGT121" "VRC07-523-LS_VRC26.25" \
        "VRC07-523-LS_PGDM1400" "VRC07-523-LS_10-1074" \
        "VRC07-523-LS_PGT121_PGDM1400" \
        "VRC01-PGDM1400-10E8v4")

outcomes=("ic80" "sens")

save_dir="/fh/fast/gilbert_p/user/bwillia2/clinical_slapnap/simulation_1b/"
mkdir -p $save_dir

io_prefix="/fh/scratch/delete90/gilbert_p/bwillia2/clinical_slapnap/output"

# Run the R script
echo -e \
    '#!/bin/bash\n Rscript simulation_1b.R ' \
    '--bnab $1 --outcome $2 --country-threshold 20 ' \
    '--output-dir $3' > run_sim1b.sh
chmod u+x run_sim1b.sh

for bnab in ${bnabs[@]}; do
    for outcome in ${outcomes[@]}; do
        echo "Running bnab: ${bnab} and outcome ${outcome}"
        full_io_prefix="${io_prefix}_${bnab}_${outcome}"
        mkdir -p $full_io_prefix
        io_file="${full_io_prefix}/slurm-%A.out"
        sbatch -A gilbert_p --mem=100G --time=1-0 -e $io_file \
            -o $io_file ./run_sim1b.sh $bnab $outcome $save_dir
    done
done
