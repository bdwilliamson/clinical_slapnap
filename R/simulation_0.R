# examine what amount of neutralization data constitutes "sufficient" to use SLAPNAP
# we'll use a plasmode simulation: subsample from VRC01, use a lasso, estimate CV-MSE

# load required functions and packages -----------------------------------------
library("here")
library("tidyverse")
library("data.table")
library("argparse")
library("SuperLearner")
library("glmnet")
library("ranger")
library("xgboost")
library("parallel")

source(here("R", "00_sim_utils.R"))
source(here("R", "00_utils.R"))

bnab <- "vrc01"
date <- "03Dec2021"
outcome_type <- "continuous"

# set up the dataset -----------------------------------------------------------
# read in the data
dat <- read_csv(here::here("dat", "slapnap_w_country_counts", paste0("slapnap_", toupper(bnab), "_", date, "_wcountry.csv")))

# redefine the outcome of interest, remove unneccessary variables
refined_dat <- dat %>%
  select(-seq.id.lanl, -seq.id.catnap, -num.seqs.in.country,
         -contains("ic50"), -contains("censored"), -country.iso) %>%
  rename(ic80 = contains("ic80")) %>%
  mutate(sens = as.numeric(ic80 < 1), .after = "ic80") %>%
  mutate(ic80 = log10(ic80)) %>% 
  select(-sens)
complete_dat <- refined_dat[complete.cases(refined_dat), ]

ns <- seq(from = 10, to = 200, by = 10)
nsim <- 2500
num_cores <- parallel::detectCores()

output_mat <- matrix(NA, nrow = length(ns) * nsim, ncol = 3)
set.seed(20230407)
seeds <- round(runif(n = length(ns), 1e3, 1e4))
cl <- parallel::makePSOCKcluster(num_cores)
parallel::clusterExport(cl = cl, varlist = ls())
parallel::clusterEvalQ(cl = cl, library("glmnet"))
parallel::clusterEvalQ(cl = cl, library("dplyr"))
for (i in 1:length(ns)) {
  clusterSetRNGStream(cl = cl, iseed = seeds[i])
  parallel::clusterExport(cl = cl, varlist = "i")
  this_output_lst <- parallel::parLapply(cl = cl, as.list(1:nsim), function(j) {
    # generate a dataset
    this_dat <- subsample_dataset(data = complete_dat, n = ns[i])
    # compute CV-MSE of a lasso predicting the outcome
    this_cvmse <- get_lasso_cvmse(data = this_dat, K = 5)
    c(ns[i], j, this_cvmse)
  })
  this_output_mat <- do.call(rbind, this_output_lst)
  output_mat[1:nsim + (i - 1) * nsim, ] <- this_output_mat
}
parallel::stopCluster(cl)
colnames(output_mat) <- c("n", "mc_id", "mse")
output_tib <- tibble::as_tibble(output_mat)
saveRDS(output_tib, file = here::here("R_output", "sim_0_results.rds"))
