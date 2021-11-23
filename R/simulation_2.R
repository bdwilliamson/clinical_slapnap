#!/usr/local/bin/Rscript
# clinical slapnap simulation 2: improving power using slapnap -----------------
# load required functions and packages -----------------------------------------
library("here")
library("tidyverse")
library("data.table")
library("argparse")
library("sievePH") # for Lunn and McNeil (1995) method

source(here("R", "00_sim_utils.R"))

# get command-line args --------------------------------------------------------
# these specify: the number of total replicates, the number of replicates per job
parser <- ArgumentParser()
parser$add_argument("--nreps-total", default = 1000, type = "double", help = "the total number of replicates")
parser$add_argument("--nreps-per-job", default = 1000, type = "double", help = "the total number of replicates per job")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1"), help = "the output directory")
args <- parser$parse_args()

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  job_id <- 1
} else {
  job_id <-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

# define static simulation parameters ------------------------------------------
nreps_per_combo <- args$nreps_total / args$nreps_per_job

# define "important" AA positions in gp120 (from AMP sieve analysis plan)
positions <- c(60, 142, 144, 147, 156, 170, 229, 230, 234, 279, 280, 317, 365, 429, 456, 458, 459, 471, 616, 824)
J <- 856

# read in the gamma values
gamma <- rep(.2, J)

# sample size
ns <- c(1e3, 5e3, 10e3, 30e3)

# set up the other simulation parameters
analyses <- c("site-scanning", "priority")
param_grid <- expand.grid(mc_id = seq_len(nreps_per_combo), n = ns, analysis = analyses)
current_dynamic_args <- param_grid[job_id, ]
print(paste0("Running n = ", current_dynamic_args$n, " for the ", current_dynamic_args$analysis, " analysis"))

# run the simulation -----------------------------------------------------------
current_seed <- job_id + current_dynamic_args$n + as.numeric(current_dynamic_args$analysis == "site-scanning") * 1e4
print(current_seed)
set.seed(current_seed)
output_lst <- lapply(as.list(1:args$nreps_per_job), function(i) {
  run_sim2_once(mc_id = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
                n = current_dynamic_args$n, J = J, gamma = gamma, delta = 1 - log(.29),
                lambda = -log(.85)/3, eos = 3, site_scanning = (current_dynamic_args$analysis == "site-scanning"),
                positions = positions)
})
output <- tibble::as_tibble(data.table::rbindlist(output_lst))
saveRDS(output, file = paste0(args$output_dir, "/output_", job_id, ".rds"))