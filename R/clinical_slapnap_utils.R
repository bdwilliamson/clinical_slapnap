# useful functions

# get the learner closest to the current date
get_learner <- function(filepath, date, outcome = "sens") {
  # list files that have "learner" but not "cv"
  all_files <- list.files(filepath, pattern = "^(?:[^c]+|c(?:$|[^v]))*$")
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