# Simulation 1b:
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
source(here("R", "03_super_learner_libraries.R"))

# get command-line args --------------------------------------------------------
# these specify the bnAb of interest, min number of readouts to specify the country(ies) of interest,
# and the outcome type (binary or continuous)
parser <- ArgumentParser()
parser$add_argument("--bnab", default = "VRC07-523-LS_PGT121", help = "the bnAb of interest")
parser$add_argument("--country-threshold", type = "double", default = "20", help = "min number of neutralization readouts for consideration")
parser$add_argument("--outcome", default = "sens", help = "the outcome of interest")
parser$add_argument("--output-dir", default = here::here("R_output", "simulation_1b"), help = "the output directory")
args <- parser$parse_args()

print(paste0("Running bnAb ", args$bnab, "; outcome ", args$outcome))

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

# get predictions from SLAPNAP--------------------------------------------------
# read in the correct model
slapnap_folder <- here::here("docker_output",
                             paste0(switch(as.numeric(outcome_type == "continuous") + 1,
                                           NULL, "continuous/"),
                                    tolower(args$bnab)))
all_learner_files <- list.files(slapnap_folder,
                           pattern = paste0("learner_", outcome_txt, "_", args$bnab))
learner_files <- all_learner_files[!grepl("cv", all_learner_files)]
dates <- as.Date(gsub(paste0("learner_", outcome_txt, "_", args$bnab, "_"), "",
                      gsub(".rds", "", learner_files)),
                    format = "%d%b%Y")
current_date <- Sys.Date()
closest_date <- which.min(current_date - dates)
slapnap_mod <- readRDS(paste0(slapnap_folder, "/", learner_files[closest_date]))

# create the correct learner/screen combinations
if (any(grepl("var_thresh", slapnap_mod$libraryNames))) {
    thresh_init <- gsub(".*?var_thresh_(.*?)", "", slapnap_mod$libraryNames)
    unique_thresh <- unique(as.numeric(thresh_init[!grepl("SL", thresh_init)]))
    for (i in seq_len(length(unique_thresh))) {
        make_screen_wrapper(var_thresh = unique_thresh[i])
    }
}
# if there were any xgboosts fit, need to "complete" them
is_xgboost <- grepl("xgboost", slapnap_mod$libraryNames)
if (any(is_xgboost)) {
    all_xgboosts <- slapnap_mod$libraryNames[is_xgboost]
    for (i in seq_len(length(all_xgboosts))) {
        slapnap_mod$fitLibrary[is_xgboost][[i]]$object <- xgb.Booster.complete(slapnap_mod$fitLibrary[is_xgboost][[i]]$object)
    }
}

# run the analyses -------------------------------------------------------------
set.seed(123125)
output <- NULL
for (i in seq_len(length(countries_of_interest))) {
  country <- countries_of_interest[i]
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
  var_theta <- est_var(dat = analysis_dataset, preds = preds, lambda = lambda_n, theta = theta)
  theta_aug <- est_theta(dat = analysis_dataset, preds = g_n, lambda = lambda_n, augmented = TRUE)
  var_theta_aug <- est_var(dat = analysis_dataset, preds = preds, lambda = lambda_n, theta = theta_aug)


  output <- bind_rows(
    output,
    tibble::tibble(bnab = args$bnab, outcome = args$outcome,
                   country = country, n_catnap = n_catnap, n_overall = n_overall,
                   augmented = c(FALSE, TRUE),
                   est = c(theta, theta_aug),
                   est_var = c(var_theta, var_theta_aug))
  )
}

saveRDS(output, paste0(args$output_dir, "/", args$bnab, "_", args$outcome, ".rds"))

print("Analysis complete.")
