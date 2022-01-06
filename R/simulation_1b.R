# Simulation 1a:
# Data analysis comparing
# (a) using sequences with data only in CATNAP to
# (b) using sequences with data in LANL, with SLAPNAP predictions

# load required functions and packages -----------------------------------------
library("here")
library("tidyverse")
library("data.table")
library("argparse")

source(here("R", "00_sim_utils.R"))

# get command-line args --------------------------------------------------------
# these specify the bnAb of interest, min number of readouts to specify the country(ies) of interest,
# and the outcome type (binary or continuous)
parser <- ArgumentParser()
parser$add_argument("--bnab", default = "PGT121", help = "the bnAb of interest")
parser$add_argument("--country_threshold", type = "double", default = "20", help = "min number of neutralization readouts for consideration")
parser$add_argument("--outcome", default = "ic80", help = "the outcome of interest")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1b"), help = "the output directory")
args <- parser$parse_args()

print(paste0("Running bnAb ", args$bnab, " for country ", args$country))

outcome_type <- ifelse(args$outcome == "ic80", "continuous", "binary")

# set up the dataset -----------------------------------------------------------
# read in the data
dat <- read_csv(here::here("data", paste0("slapnap_lanl_", args$bnab, "_wcountry.csv")))

# redefine the outcome of interest, remove unneccessary variables
refined_dat <- dat %>%
  select(-seq.id.lanl, -seq.id.catnap, -num.seqs.in.country.catnap, -num.seqs.in.country.all,
         -contains("ic50"), -contains("censored")) %>%
  rename(ic80 = !!(paste0(args$bnab, ".ic80.imputed")))

if (args$outcome == "sens") {
  refined_dat <- refined_dat %>%
    mutate(sens = ic80 < 1) %>%
    select(-ic80)
}
# define the countries of interest
all_countries <- as.list(unique(dat$country.iso))
countries_above_threshold <- unlist(lapply(all_countries, function(country) {
  this_dat <- refined_dat %>% filter(country.iso == country)
  sum(!is.na(this_dat$ic80)) > args$country_threshold
}))
countries_of_interest <- unlist(all_countries)[countries_above_threshold]

# get predictions from SLAPNAP--------------------------------------------------
# read in the correct model
all_learner_files <- list.files(here::here("docker_output",
                                      paste0(switch(as.numeric(outcome_type == "continuous") + 1, NULL, "continuous"), tolower(args$bnab))),
                           pattern = paste0("learner_", args$outcome, "_", args$bnab))
learner_files <- all_learner_files[!grepl("cv", all_learner_files)]
dates <- as.Date(gsub(paste0("learner_", args$outcome, "_", args$bnab, "_"), "", gsub(".rds", "", learner_files)),
                    format = "%d%b%Y")
current_date <- Sys.Date()
closest_date <- which.min(current_date - dates)
slapnap_mod <- readRDS(here::here("docker_output",
                                  switch(as.numeric(outcome_type == "continuous") + 1, NULL, "continuous")),
                       tolower(args$bnab),
                       learner_files[closest_date])

# run the analyses -------------------------------------------------------------
output <- NULL
for (country in countries_of_interest) {
  country_data <- refined_dat %>%
    filter(country.iso == country) %>%
    select(-country.iso)
  X <- country_data %>%
    select(-!!args$outcome)
  y <- country_data %>% pull(!!args$outcome)
  n_catnap <- sum(!is.na(y))
  n_overall <- nrow(country_data)

  # get predictions from SLAPNAP
  # preds <- predict(slapnap_mod, newdata = X, onlySL = TRUE)$pred
  preds <- predict(slapnap_mod, newdata = X)
  if (args$outcome == "ic80") {
    w <- preds
  } else {
    w <- expit(preds)
  }
  # set up the dataset
  analysis_dataset <- tibble::tibble(w = w, r = as.numeric(!is.na(y)),
                                     y = y)
  # estimate lambda
  lambda_n <- n_catnap / n_overall
  # estimate g
  g_n <- est_g(dat = analysis_dataset, type = outcome_type)
  # estimate parameter of interest
  theta <- est_theta(dat = analysis_dataset, preds = g_n, lambda = lambda_n, augmented = FALSE)
  theta_aug <- est_theta(dat = analysis_dataset, preds = g_n, lambda = lambda_n, augmented = TRUE)

  output <- bind_rows(
    output,
    tibble::tibble(bnab = args$bnab, outcome = args$outcome,
                   country = country, n_catnap = n_catnap, n_overall = n_overall,
                   augmented = c(FALSE, TRUE),
                   est = c(theta, theta_aug))
  )
}

saveRDS(output, paste0(args$output_dir, "/", args$bnab, "_", args$outcome, ".rds"))
