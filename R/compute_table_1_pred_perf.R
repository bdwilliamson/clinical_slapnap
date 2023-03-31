# compute prediction performance for each bnAb regimen (for Table 1)
# prediction performance measures: 
# CV-R^2 (continuous outcomes)
# CV-AUC, CV-Accuracy, CV-Sens/Spec/PPV/NPV @ threshold, CV-MCC, CV-Brier (binary outcomes)

# set up libraries and list of nabs for CATNAP queries
library("here")
library("tidyverse")
library("SuperLearner")
library("vimp")
source(here("R", "clinical_slapnap_utils.R"))


# load all learner objects -----------------------------------------------------
regimens <- c("vrc01", "pgt121", "vrc01-pgdm1400-10e8v4", "vrc07-523-ls", "vrc26.25",
              "vrc07-523-ls_10-1074", "vrc07-523-ls_pgdm1400", "vrc07-523-ls_pgt121", 
              "vrc07-523-ls_pgt121_pgdm1400")
folder_list_binary <- lapply(as.list(regimens), function(l) {
  here("docker_output", l)
})
folder_list_continuous <- lapply(as.list(regimens), function(l) {
  here("docker_output", "continuous", l)
})
outcomes_list_binary <- c(as.list(rep("sens", 5)), rep(list(c("estsens", "multsens")), 4))
learners_binary <- sapply(seq_len(length(folder_list_binary)), function(i) {
  get_learner(folder_list_binary[[i]], Sys.Date(), outcome = outcomes_list_binary[[i]],
              cv = TRUE)
}, simplify = FALSE)
# get prediction performance for binary outcomes -------------------------------
bin_perf <- tibble::tibble(data.table::rbindlist(lapply(1:length(learners_binary), function(l) {
  if (length(learners_binary[[l]]) == 2) {
    data.table::rbindlist(lapply(1:2, function(j) {
      get_performance_metrics(learner = learners_binary[[l]][[j]], outcome_type = "binary") %>% 
        mutate(outcome = c("estsens", "multsens")[j], regimen = regimens[l])  
    }))
  } else {
    get_performance_metrics(learner = learners_binary[[l]], outcome_type = "binary") %>% 
      mutate(outcome = "sens", regimen = regimens[l])
  }
})))

# rearrange
bin_perf <- bin_perf %>% 
  select(regimen, outcome, everything())
# look at performance at the quantile that maximized F1 score
bin_perf %>% 
  group_by(outcome, regimen) %>% 
  filter(abs(F1 - max(F1, na.rm = TRUE)) < 0.0005) %>% 
  readr::write_csv(file = here::here("R_output", "supp_table_max_f1.csv"))

# look at performance at the quantile that maximizes MCC
bin_perf %>% 
  group_by(outcome, regimen) %>% 
  filter(abs(MCC - max(MCC, na.rm = TRUE)) < 0.0005) %>% 
  readr::write_csv(file = here::here("R_output", "supp_table_max_mcc.csv"))

# print out a table of binary performance metrics, for each regimen and outcome
for (i in seq_len(length(regimens))) {
  bin_perf %>% 
    filter(regimen == regimens[i]) %>% 
    readr::write_csv(file = here::here("R_output", paste0("supp_table_", regimens[i], ".csv")))
}
