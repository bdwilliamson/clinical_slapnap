#!/bin/bash

# all of the nabs in Table 1
ALL_NABS=("VRC01" "VRC07-523-LS" "PGT121" \
          "VRC07-523-LS+PGDM1400" \
          "VRC07-523-LS+10-1074" \
          "VRC07-523-LS+PGT121" \
          "VRC07-523-LS+PGT121+PGDM1400" \
          "VRC01/PGDM1400-10E8v4")
ALL_NABS=("VRC07-523-LS+PGT121" \
          "VRC07-523-LS+PGT121+PGDM1400")
# Notes:
# 10E8v4/iMab bi-specific is not in CATNAP yet; NCT03875209
# PGT121.414.LS also not in CATNAP; NCT04212091
# CAP256V2LS also not in CATNAP

# for each nab, run slapnap
for NAB in ${ALL_NABS[@]}
do
  SLAPNAP_NAB=$(echo $NAB | sed 's/+/;/g' | sed 's/ //g')
  DIR_NAB=$(echo $NAB | tr '[:upper:]' '[:lower:]' | sed 's/+/_/g' | sed 's/ //g')
  mkdir -p ./docker_output/$DIR_NAB
  if grep -q ";" <<< "$SLAPNAP_NAB"; then
      outcomes="estsens;multsens"
  else
      outcomes="sens"
  fi
  # run variable importance for a subset of mAbs
  if grep -q "VRC07-523-LS;PGT121" <<< "$SLAPNAP_NAB"; then
      importance_grp="cond;marg"
      importance_ind="marg;pred"
  else
      importance_grp=""
      importance_ind=""
  fi
  sudo docker run \
    -v "$(pwd)"/docker_output/$DIR_NAB/:/home/output \
    -e learners="rf;lasso;xgboost" \
    -e cvperf="TRUE" \
    -e cvtune="TRUE" \
    -e nab=$SLAPNAP_NAB \
    -e outcomes=$outcomes \
    -e sens_thresh="1" \
    -e binary_outcomes="ic80" \
    -e importance_grp=$importance_grp \
    -e importance_ind=$importance_ind \
    -e var_thresh="0;4" \
    -e return="report;data;learner;figures;vimp" \
    -e nfolds="5" \
    slapnap
done
