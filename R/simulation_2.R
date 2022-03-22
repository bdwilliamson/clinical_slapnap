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
parser$add_argument("--nreps-total", default = 1000, type = "double", help = "the total number of replicates")
parser$add_argument("--nreps-per-job", default = 1, type = "double", help = "the total number of replicates per job")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1"), help = "the output directory")
args <- parser$parse_args()

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  job_id <- 3001
} else {
  job_id <-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

# define static simulation parameters ------------------------------------------
nreps_per_combo <- args$nreps_total / args$nreps_per_job

# read in data from AMP
sens_data_from_amp <- readr::read_csv(here::here("dat", "resis_resid_full_v2.csv"))

# define "important" AA positions in gp120 (from AMP sieve analysis plan)
positions <- sens_data_from_amp %>% filter(preselected.site == 1) %>% pull(pos)
# from Craig: the positions passing minimum variability filter
minvar_positions <- sens_data_from_amp$pos
all_positions <- sens_data_from_amp %>% filter(pos %in% minvar_positions) %>% pull(pos)
J <- nrow(sens_data_from_amp)

# proportions of infected subjects in amp with putative sensitive genotype at each position, each tx group
gammas <- sens_data_from_amp %>%
  select(pos, prop.pooled.placebo.w.sens) %>%
  rename(prop_placebo = prop.pooled.placebo.w.sens) %>%
  filter(pos %in% minvar_positions)

# note that at site 230, residue is D, which was actually resistant
gammas[gammas$pos == 230, ]$prop_placebo <- 1 - gammas[gammas$pos == 230, ]$prop_placebo

# assume PE(overall) = 0.7
pe_overall <- 0.7

# prevention efficacy given a sensitive residue at the given position;
# assume PE(sensitive at position j) = PE(other at position j) = 0.7 at all positions besides 230
# assume PE(resistant at 230) = 0
pe_230 <- 1 - exp(log(1 - pe_overall) / gammas$prop_placebo[gammas$pos == 230])
pe_sites <- rep(0.7, nrow(gammas))
pe_sites[gammas$pos == 230] <- pe_230

# rate of loss-to-followup per year; assume 10% (observed 9.4% in HVTN 704, 6.3% in 703)
# q <- 769 / 4611
q <- 0.1

# sample size; from Gilbert (2019), Table 2
ns <- c(2071, 4141, 12422)
# incidence among placebo-arm subjects; rate per person-year
lambda_0s <- c(.037, .018, .006)

# set up the other simulation parameters
analyses <- c("site-scanning", "priority")
param_grid <- expand.grid(mc_id = seq_len(nreps_per_combo), n = ns, analysis = analyses)
param_grid$lambda_0 <- case_when(
  param_grid$n == ns[1] ~ lambda_0s[1],
  param_grid$n == ns[2] ~ lambda_0s[2],
  param_grid$n == ns[3] ~ lambda_0s[3]
)
current_dynamic_args <- param_grid[job_id, ]
print(paste0("Running n = ", current_dynamic_args$n, " for the ", current_dynamic_args$analysis, " analysis"))

# run the simulation -----------------------------------------------------------
current_seed <- job_id + current_dynamic_args$n + as.numeric(current_dynamic_args$analysis == "site-scanning") * 1e4
print(current_seed)
set.seed(current_seed)
# approx. 38 seconds per replication for site-scanning (n = 4000).
# approx. 1 second per replication for priority (n = 4000).
system.time(
  output_lst <- lapply(as.list(1:args$nreps_per_job), function(i) {
    run_sim2_once(mc_id = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
                  n = current_dynamic_args$n, all_positions = all_positions, gamma = gammas,
                  lambda_0 = current_dynamic_args$lambda_0, beta = log(1 - pe_overall),
                  alpha = log(1 - pe_sites) - log(1 - 0), 
                  eos = 365 * 2, q = q,
                  site_scanning = (current_dynamic_args$analysis == "site-scanning"),
                  positions = positions, position = which(gammas$pos == 230),
                  minvar_screen = 0)
  })
)
output <- tibble::as_tibble(data.table::rbindlist(output_lst))
# for debugging: 
# round(colMeans(output))
output %>%
  group_by(position) %>%
  summarize(power = mean(reject, na.rm = TRUE))
saveRDS(output, file = paste0(args$output_dir, "/output_", job_id, ".rds"))
