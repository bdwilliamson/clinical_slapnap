#!/bin/bash

# run the entire simulation 1
ml fhR/4.0.2-foss-2019b
ml jbigkit

# all bnabs
bnabs=("VRC01" "VRC07-523-LS" "PGT121" "VRC26.25" "PGDM1400" \
        "VRC07-523-LS_PGT121" "VRC07-523-LS_VRC26.25" \
        "VRC07-523-LS_PGDM1400" "VRC07-523-LS_10-1074" \
        "VRC07-523-LS_PGT121_PGDM1400" \
        "VRC01-PGDM1400-10E8v4")
outcomes=("ic80" "sens")
simplifiers=(0 1)

nreps_per_job=1000
save_dir="/fh/fast/gilbert_p/user/bwillia2/clinical_slapnap/simulation_1/"
mkdir -p $save_dir

for bnab in ${bnabs[@]}; do
    for outcome in ${outcomes[@]}; do
        for simplify in ${simplifiers[@]}; do
            io_prefix="clinical_slapnap/output_${bnab}_${outcome}"
            ./run_simulation_1.sh $bnab $outcome 1000 $nreps_per_job $save_dir $io_prefix $simplify
        done
    done
done
