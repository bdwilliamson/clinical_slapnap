# Simulation 1a:
# Data analysis comparing
# (a) using sequences with data only in CATNAP to
# (b) using sequences with data in LANL, with SLAPNAP predictions

# load required functions and packages -----------------------------------------
library("here")
library("tidyverse")
library("data.table")
library("argparse")
library("SuperLearner")
library("glmnet")
library("ranger")
library("xgboost")

source(here("R", "00_sim_utils.R"))
source(here("R", "00_utils.R"))

# get command-line args --------------------------------------------------------
# these specify the bnAb of interest, min number of readouts to specify the country(ies) of interest,
# and the outcome type (binary or continuous)
parser <- ArgumentParser()
parser$add_argument("--bnab", default = "PGT121", help = "the bnAb of interest")
parser$add_argument("--country-threshold", type = "double", default = "20", help = "min number of neutralization readouts for consideration")
parser$add_argument("--outcome", default = "ic80", help = "the outcome of interest")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1b"), help = "the output directory")
args <- parser$parse_args()

print(paste0("Running bnAb ", args$bnab))

outcome_type <- ifelse(args$outcome == "ic80", "continuous", "binary")

# set up the dataset -----------------------------------------------------------
# read in the data
dat <- read_csv(here::here("dat", paste0("slapnap_lanl_", args$bnab, "_wcountry.csv")))

# redefine the outcome of interest, remove unneccessary variables
if (!any(grepl("_", args$bnab))) {
  # handle single bnAb
  refined_dat <- dat %>%
    select(-seq.id.lanl, -seq.id.catnap, -num.seqs.in.country.catnap, -num.seqs.in.country.all,
           -seq.in.catnap, -seq.in.lanl.2019, -contains("ic50"), -contains("censored")) %>%
    rename(ic80 = contains("ic80")) %>%
    mutate(sens = as.numeric(ic80 < 1), .after = "ic80")
  n_ab <- 1
} else {
  # handle bnAb regimen
  refined_dat <- dat %>%
    select(-seq.id.lanl, -seq.id.catnap, -num.seqs.in.country.catnap, -num.seqs.in.country.all,
           -seq.in.catnap, -seq.in.lanl.2019, -contains("ic50"), -contains("ic80"),
           -contains("iip"), -contains("dichotomous")) %>%
    mutate(ic80 = dat %>% pull(pc.ic80),
           sens = as.numeric(ic80 < 1), .after = "country.iso")
  n_ab <- 2
}

# define the countries of interest
all_countries <- as.list(unique(dat$country.iso))
countries_above_threshold <- unlist(lapply(all_countries, function(country) {
  this_dat <- refined_dat %>% filter(country.iso == country)
  sum(!is.na(this_dat$ic80)) > args$country_threshold
}))
countries_of_interest <- unlist(all_countries)[countries_above_threshold]

# remove ic80 if outcome is sens; otherwise, remove sens
if (args$outcome == "sens") {
  refined_dat <- refined_dat %>%
    select(-ic80)
} else {
  refined_dat <- refined_dat %>%
    select(-sens)
}

if (n_ab == 1) {
  outcome_txt <- args$outcome
} else {
  if (args$outcome == "sens") {
    outcome_txt <- "estsens"
  } else {
    outcome_txt <- "ic80"
  }
}
b <- 5000

# get predictions from SLAPNAP--------------------------------------------------
# read in the correct model
slapnap_folder <- here::here("docker_output",
                             paste0(switch(as.numeric(outcome_type == "continuous") + 1, NULL, "continuous"), "/", tolower(args$bnab)))
all_learner_files <- list.files(slapnap_folder,
                           pattern = paste0("learner_", args$outcome, "_", args$bnab))
learner_files <- all_learner_files[!grepl("cv", all_learner_files)]
dates <- as.Date(gsub(paste0("learner_", args$outcome, "_", args$bnab, "_"), "", gsub(".rds", "", learner_files)),
                    format = "%d%b%Y")
current_date <- Sys.Date()
closest_date <- which.min(current_date - dates)
slapnap_mod <- readRDS(paste0(slapnap_folder, "/", learner_files[closest_date]))

# estimate E(Y | W) based on all of slapnap preds
all_preds <- slapnap_mod$SL.predict
dat_for_g <- tibble::tibble(w = all_preds, y = slapnap_mod$Y) %>%
  mutate(r = as.numeric(!is.na(y)))
g_mod <- glm(y ~ w, data = dat_for_g,
             family = switch(as.numeric(outcome_type == "binary") + 1, gaussian(),
                             binomial()))


# run the analyses -------------------------------------------------------------
set.seed(4747)
seeds <- round(runif(length(countries_of_interest), 1e4, 1e5))
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
  preds <- predict(slapnap_mod, newdata = X, onlySL = TRUE)$pred
  # if (args$outcome == "ic80") {
  #   w <- preds
  # } else {
  #   w <- expit(preds)
  # }
  w <- preds
  # set up the dataset
  analysis_dataset <- tibble::tibble(w = w, r = as.numeric(!is.na(y)),
                                     y = y)
  set.seed(seeds[country == countries_of_interest])
  # estimate lambda
  lambda_n <- n_catnap / n_overall
  # estimate g
  g_n <- predict(g_mod, newdata = analysis_dataset, type = "response")
  # g_n <- est_g(dat = analysis_dataset, type = outcome_type)
  # estimate parameter of interest
  theta <- est_theta(dat = analysis_dataset, preds = g_n, lambda = lambda_n, augmented = FALSE)
  theta_aug <- est_theta(dat = analysis_dataset, preds = g_n, lambda = lambda_n, augmented = TRUE)
  boot_theta <- boot::boot(data = analysis_dataset, statistic = sim1b_boot_stat,
                           R = b, sim = "ordinary", stype = "i",
                           augmented = FALSE, g_n = g_mod,
                           outcome_type = outcome_type)
  boot_ci_theta <- boot::boot.ci(boot_theta, conf = 0.95, type = "perc")
  boot_theta_aug <- boot::boot(data = analysis_dataset, statistic = sim1b_boot_stat,
                               R = b, sim = "ordinary", stype = "i",
                               augmented = TRUE, g_n = g_mod,
                               outcome_type = outcome_type)
  boot_ci_theta_aug <- boot::boot.ci(boot_theta_aug, conf = 0.95, type = "perc")
  if (is.null(boot_ci_theta)) {
    boot_ci_theta <- list("percent" = matrix(c(0, 0, 0, theta, theta), nrow = 1))
  }
  if (is.null(boot_ci_theta_aug)) {
    boot_ci_theta_aug <- list("percent" = matrix(c(0, 0, 0, theta_aug, theta_aug), nrow = 1))
  }
  cis <- rbind(boot_ci_theta$percent[, 4:5], boot_ci_theta_aug$percent[, 4:5])

  output <- bind_rows(
    output,
    tibble::tibble(bnab = args$bnab, outcome = args$outcome,
                   country = country, n_catnap = n_catnap, n_overall = n_overall,
                   augmented = c(FALSE, TRUE),
                   est = c(theta, theta_aug),
                   ci_ll = cis[, 1], ci_ul = cis[, 2]) %>%
      mutate(width = ci_ul - ci_ll,
             sq_width = width ^ 2)
  )
}

saveRDS(output, paste0(args$output_dir, "/", args$bnab, "_", args$outcome, ".rds"))
