#!/usr/local/bin/Rscript
# clinical slapnap simulation 2: improving power using slapnap -----------------
# load required functions and packages -----------------------------------------
library("here")
library("tidyverse")
library("data.table")
library("argparse")
library("sievePH") # for Michal's method
library("survival") # for Lunn & McNeil (1995) method

source(here("R", "00_sim_utils.R"))
source(here("R", "lunnMcneil.R"))

# get command-line args --------------------------------------------------------
# these specify: the number of total replicates, the number of replicates per job
parser <- ArgumentParser()
parser$add_argument("--analysis", default = "priority", help = "the analysis to run")
parser$add_argument("--nreps-total", default = 1000, type = "double", help = "the total number of replicates")
parser$add_argument("--nreps-per-job", default = 500, type = "double", help = "the total number of replicates per job")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1"), help = "the output directory")
args <- parser$parse_args()

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  job_id <- 1 # n = 2071, priority
  # job_id <- 4001 # n = 4141, priority
  # job_id <- 5001 # n = 12422, priority
} else {
  job_id <-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

# define static simulation parameters ------------------------------------------
nreps_per_combo <- args$nreps_total / args$nreps_per_job

# read in data from AMP
sens_data_from_amp <- readr::read_csv(here::here("dat", "resis_resid_full_v3.csv"))

# define "important" AA positions in gp120 (from AMP sieve analysis plan)
positions <- sens_data_from_amp %>% filter(preselected.site == 1) %>% pull(pos)
# the positions passing minimum variability filter (>= 4 sequences that differ from most common
# at a given site)
sens_data_from_amp_minvar <- sens_data_from_amp %>%
  filter(minvar.pooled >= 4)
all_positions <- sens_data_from_amp_minvar %>% pull(pos)
J <- nrow(sens_data_from_amp_minvar)

# proportions of infected subjects in amp with putative sensitive genotype at each position, each tx group
gammas <- sens_data_from_amp_minvar %>%
  select(pos, prop.pooled.placebo.w.sens) %>%
  rename(prop_placebo = prop.pooled.placebo.w.sens)

# note that at site 230, residue is D, which was actually resistant
gammas[gammas$pos == 230, ]$prop_placebo <- 1 - gammas[gammas$pos == 230, ]$prop_placebo

# assume PE(overall) = 0.7
pe_overall <- 0.7

# rate of loss-to-followup per year; assume 10% (observed 9.4% in HVTN 704, 6.3% in 703)
q <- 0.1

# sample size; from Gilbert (2019), Table 2
ns <- c(2071, 4141, 12422)
# incidence among placebo-arm subjects; rate per person-year
lambda_0s <- c(.03, .018, .006) # from Gilbert (2019), but we need smaller for higher PE

# set up the other simulation parameters
param_grid <- expand.grid(mc_id = seq_len(nreps_per_combo), n = ns,
                          pe_0 = seq(0, .7, .7 / 3))
param_grid$lambda_0 <- case_when(
  param_grid$n == ns[1] ~ lambda_0s[1],
  param_grid$n == ns[2] ~ lambda_0s[2],
  param_grid$n == ns[3] ~ lambda_0s[3]
)
# prevention efficacy given a sensitive residue at the given position;
param_grid$pe_1 <- get_pe_1(pe_overall = pe_overall, pe_0 = param_grid$pe_0,
                            gamma = gammas$prop_placebo[gammas$pos == 230])
current_dynamic_args <- param_grid[job_id, ]
print(paste0("Running n = ", current_dynamic_args$n, " for the ", args$analysis, " analysis, pe_0 = ", current_dynamic_args$pe_0))

# run the simulation -----------------------------------------------------------
current_seed <- job_id + current_dynamic_args$n + as.numeric(args$analysis == "site-scanning") * 1e4
print(current_seed)
set.seed(current_seed)
if (.Platform$OS.type == "windows") {
  library("parallel")
  num_cores <- parallel::detectCores()
  cl <- parallel::makePSOCKcluster(num_cores)
  parallel::clusterExport(cl = cl, varlist = ls())
  parallel::clusterEvalQ(cl = cl, library("survival"))
  parallel::clusterEvalQ(cl = cl, library("sievePH"))
  parallel::clusterEvalQ(cl = cl, library("dplyr"))
  clusterSetRNGStream(cl = cl, iseed = current_seed)
  start <- Sys.time()
  output_lst <- parallel::parLapply(cl = cl,
    as.list(1:args$nreps_per_job), function(i) {
      run_sim2_once(mc_id = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
                    n = current_dynamic_args$n, all_positions = all_positions, gamma = gammas,
                    lambda_0 = current_dynamic_args$lambda_0, beta = log(1 - pe_overall),
                    pe_1 = current_dynamic_args$pe_1, pe_0 = current_dynamic_args$pe_0,
                    eos = 365 * 2, q = q,
                    site_scanning = (args$analysis == "site-scanning"),
                    positions = positions, position = which(gammas$pos == 230),
                    minvar_screen = 0, package = "none")
    }
  )
  # output_lst <- parallel::parLapply(
  #   cl = cl, as.list(1:args$nreps_per_job), function(i) {
  #     run_sim2_simple_once(mc_id = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
  #                          n = current_dynamic_args$n, lambda_0 = current_dynamic_args$lambda_0,
  #                          alpha = log(1 - .7), eos = 365 * 2, q = q)
  #   }
  # )
  end <- Sys.time()
  cat("Elapsed time: ", format(end - start), "\n")
  parallel::stopCluster(cl)
} else {
  start <- Sys.time()
  output_lst <- lapply(as.list(1:args$nreps_per_job), function(i) {
    run_sim2_once(mc_id = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
                  n = current_dynamic_args$n, all_positions = all_positions, gamma = gammas,
                  lambda_0 = current_dynamic_args$lambda_0, beta = log(1 - pe_overall),
                  pe_1 = current_dynamic_args$pe_1, pe_0 = current_dynamic_args$pe_0,
                  eos = 365 * 2, q = q,
                  site_scanning = (args$analysis == "site-scanning"),
                  positions = positions, position = which(gammas$pos == 230),
                  minvar_screen = 0, package = "none")
  })
  end <- Sys.time()
  cat("Elapsed time: ", format(end - start), "\n")
}
output <- tibble::as_tibble(data.table::rbindlist(output_lst))
# for debugging:
# round(colMeans(output))
output %>%
  group_by(position) %>%
  summarize(sieve_power = mean(reject, na.rm = TRUE),
            unadj_sieve_power = mean(unadjusted_p_val < 0.05, na.rm = TRUE),
            overall_power = mean(hr_p < 0.05, na.rm = TRUE),
            n_events = mean(n_events, na.rm = TRUE),
            n00 = mean(n00), n01 = mean(n01), n10 = mean(n10), n11 = mean(n11))
saveRDS(output, file = paste0(args$output_dir, "/output_", args$analysis, "_", job_id, ".rds"))
