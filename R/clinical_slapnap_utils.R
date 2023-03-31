# useful functions

# get the learner closest to the current date
get_learner <- function(filepath, date, outcome = "sens", cv = FALSE) {
  # list files that have "learner" but not "cv"
  if (!cv) {
    all_files <- list.files(filepath, pattern = "^(?:[^c]+|c(?:$|[^v]))*$")  
  } else {
    all_files <- list.files(filepath, pattern = "cv")
  }
  
  learner_files <- all_files[grepl("learner", all_files)]
  # do this for each unique outcome
  if (length(outcome) > 1) {
    learner <- vector("list", length = length(outcome))
    for (i in seq_len(length(outcome))) {
      these_learner_files <- learner_files[grepl(outcome[i], learner_files)]
      # find the date closest to now
      dates <- as.Date(unlist(lapply(
        strsplit(these_learner_files, "_", fixed = TRUE), function(x) gsub(".rds", "", tail(x, n = 1))
      )), format = "%d%b%Y")
      closest_date <- which.min(date - dates)
      learner[[i]] <- readRDS(paste0(filepath, "/", these_learner_files[closest_date]))
    }
  } else {
    # find the date closest to now
    dates <- as.Date(unlist(lapply(
      strsplit(learner_files, "_", fixed = TRUE), function(x) gsub(".rds", "", tail(x, n = 1))
    )), format = "%d%b%Y")
    closest_date <- which.min(date - dates)
    learner <- readRDS(paste0(filepath, "/", learner_files[closest_date]))
  }
  learner
}

get_f_score <- function(precision, recall, beta = 1) {
  (1 + beta ^ 2) * precision * recall / (beta ^ 2 * precision + recall)
}

# get performance metrics for a given learner
# @param learner a cross-validated SL object
get_performance_metrics <- function(learner, outcome_type = "continuous") {
  these_folds <- vimp::get_cv_sl_folds(cv_sl_folds = learner$folds)
  V <- length(unique(these_folds))
  these_r2s <- sapply(seq_len(length(unique(these_folds))), function(v) {
    vimp::measure_r_squared(fitted_values = learner$SL.predict[these_folds == v], 
                            y = learner$Y[these_folds == v])  
  }, simplify = FALSE)
  these_r2s <- unlist(lapply(these_r2s, function(r2) r2$point_est))
  this_r2 <- mean(these_r2s)
  if (outcome_type == "continuous") {
    perf <- tibble::tibble("perf_measure" = "R2", "perf" = this_point_est,
                           "cutoff" = NA)
  } else {
    quantiles <- seq(1, 99) / 100
    cutoff_dependent_perf <- data.frame("quantile" = quantiles,
                                        "cutoff" = NA,
                                        "Sensitivity" = NA,
                                        "Specificity" = NA,
                                        "PPV" = NA, "NPV" = NA,
                                        "F1" = NA, "F0.5" = NA,
                                        "Accuracy" = NA, "MCC" = NA)
    for (i in 1:length(quantiles)) {
      this_perf <- NULL
      for (v in seq_len(V)) {
        cutoff <- quantile(learner$SL.predict[these_folds != v], p = quantiles[i])
        if (is.null(this_perf)) {
          this_perf <- get_fold_perf(y = learner$Y[these_folds == v], 
                                     preds = learner$SL.predict[these_folds == v],
                                     cutoff = cutoff)  
        } else {
          this_perf <- dplyr::left_join(this_perf, get_fold_perf(y = learner$Y[these_folds == v], 
                                                                 preds = learner$SL.predict[these_folds == v],
                                                                 cutoff = cutoff),
                                        by = "perf_measure")
        }
      }
      cutoff_dependent_perf[i, ] <- c(quantiles[i], cutoff, 
                                      rowMeans(this_perf %>% select(-perf_measure)))
    }
    auc <- cvAUC::cvAUC(predictions = learner$SL.predict,
                        labels = learner$Y, 
                        folds = these_folds)$fold.AUC
    brier <- this_r2
    perf <- cutoff_dependent_perf
    perf$AUC <- auc
    perf$Brier <- brier
  }
  perf
}

get_fold_perf <- function(y, preds, cutoff) {
  TP <- sum(y == 1 & preds >= cutoff)
  FP <- sum(y == 0 & preds >= cutoff)
  TN <- sum(y == 0 & preds < cutoff)
  FN <- sum(y == 1 & preds < cutoff)
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  F1 <- (1 + 1 ^ 2) * TP / ( (1 + 1 ^ 2) * TP + 1 ^ 2 * FN + FP)
  F05 <- (1 + 0.5 ^ 2) * TP / ( (1 + 0.5 ^ 2) * TP + 0.5 ^ 2 * FN + FP)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  mcc <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  return(tibble::tibble("perf_measure" = c("Sensitivity", "Specificity", "PPV", "NPV", "F1", "F0.5", "Accuracy", "MCC"),
                        "perf" = c(sens, spec, PPV, NPV, F1, F05, accuracy, mcc)))
}
