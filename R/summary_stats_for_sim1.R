# create a table with the mean and variance of the outcome, 
# along with R-squared / (estimated) sensitivity for each bnab and outcome
library("tidyverse")
bnabs <- c("vrc01", "pgt121", "vrc07-523-ls", "vrc07-523-ls_10-1074",
           "vrc07-523-ls_pgt121", "vrc07-523-ls_pgdm1400", "vrc07-523-ls_pgt121_pgdm1400",
           "vrc01/pgdm1400-10e8v4")
# note that this goes (ic80, sens) for each bnab
ests <- c(.345, .744, .571, .85, .193, .728, .319, .783, .316, .768,
          .255, .638, .181, .73, .254, .81)

output <- tibble::tibble(bnab = rep(bnabs, each = 2),
                         outcome = rep(c("ic80", "sens"), length(bnabs)),
                         est = ests,
                         mn_y = rep(NA_real_, length(bnabs) * 2),
                         var_y = rep(NA_real_, length(bnabs) * 2))
for (i in 1:nrow(output)) {
  bnab <- output$bnab[i]
  this_folder <- here::here("docker_output", 
                            paste0(switch((output$outcome[i] == "ic80") + 1, NULL, "continuous/"),
                            bnab))
  this_datafile <- list.files(this_folder, pattern = ".csv")
  data <- readr::read_csv(paste0(this_folder, "/", this_datafile))
  this_outcome <- ifelse(output$outcome[i] == "ic80", "ic80", 
                         ifelse(grepl("_", bnab, fixed = TRUE), "estsens", "sens"))
  y <- data %>% pull(!!this_outcome)
  output[i, ]$mn_y <- mean(y, na.rm = TRUE)
  output[i, ]$var_y <- var(y, na.rm = TRUE)
}

saveRDS(output, file = here::here("R_output", "summary_statistics_simulation_1.rds"))
