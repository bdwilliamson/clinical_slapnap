# clinical slapnap simulation 1: relative efficiency using slapnap -------------
# load required functions and packages -----------------------------------------
library("here")
library("tidyverse")
library("data.table")
library("argparse")

source(here("R", "00_sim_utils.R"))

# get command-line args --------------------------------------------------------
# these specify: the bnAb of interest, the number of total replicates, the number of replicates per job, the type of outcome (binary or continuous)
parser <- ArgumentParser()
parser$add_argument("--bnab", default = "VRC01", help = "the bnAb of interest")
parser$add_argument("--outcome", default = "ic80", help = "the outcome of interest")
parser$add_argument("--simplify", default = 1, type = "double", help = "should we use the same mean for Y and W?")
parser$add_argument("--nreps-total", default = 1000, type = "double", help = "the total number of replicates")
parser$add_argument("--nreps-per-job", default = 1000, type = "double", help = "the total number of replicates per job")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1"), help = "the output directory")
args <- parser$parse_args()

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  job_id <- 1
} else {
  # modify the following line if you are not running on a Slurm system
  job_id <-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

# set up the varying simulation parameters
n <- 1000
epsilons <- c(0.5, 1, 2)
nreps_per_combo <- args$nreps_total / args$nreps_per_job
param_grid <- expand.grid(mc_id = 1:nreps_per_combo, n = n, epsilon = epsilons)
current_dynamic_args <- param_grid[job_id, ]
print(paste0("Running epsilon = ", current_dynamic_args$epsilon, " for bnAb regimen ", args$bnab))

# set up the fixed simulation parameters
summary_statistics <- readRDS(here::here("R_output", "summary_statistics_simulation_1.rds"))
these_summary_statistics <- summary_statistics %>%
  filter(bnab == tolower(args$bnab), outcome == args$outcome)

# set up the parameters list
if (args$outcome == "ic80") {
    params_lst <- list(mu1 = these_summary_statistics$mn_y,
                       var1 = these_summary_statistics$var_y,
                       mu2 = switch(as.numeric(args$simplify) + 1, these_summary_statistics$mn_w, these_summary_statistics$mn_y),
                       var2 = switch(as.numeric(args$simplify) + 1, these_summary_statistics$var_w, these_summary_statistics$var_y))
} else {
    params_lst <- list(mu0 = -0.32, sigma0 = 0.2, sigma1 = 0.2,
                       p_y = these_summary_statistics$mn_y)
}

# run the simulation -----------------------------------------------------------
current_seed <- job_id
print(current_seed)
set.seed(current_seed)
output_lst <- lapply(as.list(1:args$nreps_per_job), function(i) {
    get_ests(
        mc_id = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
        n = current_dynamic_args$n,
        epsilon = current_dynamic_args$epsilon,
        point_est = these_summary_statistics$est,
        datatype = ifelse(args$outcome == "ic80", "continuous", "binary"),
        params = params_lst
    )
})
output <- tibble::as_tibble(data.table::rbindlist(output_lst)) %>%
  mutate(bnab = args$bnab, simplified = args$simplify, .before = "mc_id")
saveRDS(output, file = paste0(args$output_dir, "/output_", args$outcome, "_",
                              tolower(args$bnab), "_", args$simplify, "_", job_id, ".rds"))
print("Analysis complete!")
