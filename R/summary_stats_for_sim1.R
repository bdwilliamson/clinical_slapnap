# create a table with the mean and variance of the outcome, 
# along with R-squared / (estimated) sensitivity for each bnab and outcome
library("tidyverse")
bnabs <- c("vrc01", "vrc07-523-ls", "pgt121", "vrc26.25", "pgdm1400", 
           "vrc07-523-ls_pgt121", "vrc07-523-ls_vrc26.25", "vrc07-523-ls_pgdm1400", "vrc07-523-ls_10-1074", 
           "vrc07-523-ls_pgt121_pgdm1400",
           "vrc01-pgdm1400-10e8v4")
# note that this goes (ic80, sens) for each bnab
ests <- c(.345, .744, # vrc01
          .193, .728, # vrc07-523-ls
          .571, .85, # pgt121
          .53, .867, # vrc26.25
          .501, .873, # pgdm1400
          .316, .768, # vrc07-523-ls + pgt121
          .373, .689, # vrc07-523-ls + vrc26.25
          .255, .638, # vrc07-523-ls + pgdm1400
          .319, .783, # vrc07-523-ls + 10-1074
          .181, .73, # vrc07-523-ls + pgt121 + pgdm1400
          .254, .81) # vrc01/pgdm1400-10e8v4

output <- tibble::tibble(bnab = rep(gsub("/", "-", bnabs), each = 2),
                         outcome = rep(c("ic80", "sens"), length(bnabs)),
                         est = ests,
                         mn_y = rep(NA_real_, length(bnabs) * 2),
                         var_y = rep(NA_real_, length(bnabs) * 2),
                         mn_w = rep(NA_real_, length(bnabs) * 2),
                         var_w = rep(NA_real_, length(bnabs) * 2))
for (i in 1:nrow(output)) {
  bnab <- bnabs[which(gsub("/", "-", bnabs) == output$bnab[i])]
  this_folder <- here::here("docker_output", 
                            paste0(switch((output$outcome[i] == "ic80") + 1, NULL, "continuous/"),
                            bnab))
  this_datafile <- list.files(this_folder, pattern = ".csv")[length(list.files(this_folder, pattern = ".csv"))]
  this_cvlearnerfile <- list.files(this_folder, pattern = "cvlearner")[length(list.files(this_folder, pattern = "cvlearner"))]
  data <- readr::read_csv(paste0(this_folder, "/", this_datafile))
  cvlearner <- readRDS(paste0(this_folder, "/", this_cvlearnerfile))
  this_outcome <- ifelse(output$outcome[i] == "ic80", "ic80", 
                         ifelse(grepl("_", bnab, fixed = TRUE), "estsens", "sens"))
  y <- data %>% pull(!!this_outcome)
  output[i, ]$mn_y <- mean(y, na.rm = TRUE)
  output[i, ]$var_y <- var(y, na.rm = TRUE)
  output[i, ]$mn_w <- mean(cvlearner$SL.predict)
  output[i, ]$var_w <- var(cvlearner$SL.predict)
}

saveRDS(output, file = here::here("R_output", "summary_statistics_simulation_1.rds"))
