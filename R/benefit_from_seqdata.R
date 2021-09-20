# compute the benefit (reduction in variance, CI width) from using 
# additional sequence data via SLAPNAP over just the data with measured IC80 in 
# CATNAP for each bnAb

# load required libraries ------------------------------------------------------
library("tidyverse")
library("here")
library("SuperLearner")
library("vimp")
source(here("R", "clinical_slapnap_utils.R"))

# load all learner objects -----------------------------------------------------
regimens <- c("vrc01", "pgt121", "vrc01-pgdm1400-10e8v4", "vrc07-523-ls", 
              "vrc07-523-ls_10-1074", "vrc07-523-ls_pgdm1400", "vrc07-523-ls_pgt121", 
              "vrc07-523-ls_pgt121_pgdm1400")
folder_list <- lapply(as.list(regimens), function(l) {
  here("docker_output", l)
})
outcomes_list <- c(as.list(rep("sens", 4)), rep(list(c("estsens", "multsens")), 4))
learners <- sapply(seq_len(length(folder_list)), function(i) {
  get_learner(folder_list[[i]], Sys.Date(), outcome = outcomes_list[[i]])
})
# split up into single nAbs, combo nAbs (and est, mult sens)
single_nabs <- learners[1:4]
estsens <- lapply(learners[5:8], function(l) l[[1]])
multsens <- lapply(learners[5:8], function(l) l[[2]])

# get R-squared for each -------------------------------------------------------
single_r2 <- unlist(lapply(single_nabs, function(l) vimp::measure_r_squared(y = l$Y, fitted_values = l$SL.predict)$point_est))
estsens_r2 <- unlist(lapply(estsens, function(l) vimp::measure_r_squared(y = l$Y, fitted_values = l$SL.predict)$point_est))
multsens_r2 <- unlist(lapply(multsens, function(l) vimp::measure_r_squared(y = l$Y, fitted_values = l$SL.predict)$point_est))
r2_tib <- tibble(bnab = regimens, estsens = c(single_r2, estsens_r2), multsens = c(rep(NA, 4), multsens_r2))
readr::write_csv(r2_tib, here("R_output", "table_1_r2s.csv"))
