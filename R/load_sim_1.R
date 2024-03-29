#!/usr/local/bin/Rscript
library("here")
library("tidyverse")
library("data.table")
library("cowplot")
theme_set(theme_cowplot())

summary_statistics <- readRDS(here::here("R_output", "summary_statistics_simulation_1.rds"))

# read in results of the simulation, create a table
all_output_files <- list.files(here::here("R_output", "simulation_1"))
all_output <- tibble::tibble(data.table::rbindlist(
  lapply(as.list(here::here("R_output", "simulation_1", all_output_files)), readRDS)
)) %>% 
  left_join(summary_statistics %>% 
              mutate(datatype = ifelse(outcome == "ic80", "continuous", "binary"),
                     bnab = toupper(bnab),
                     bnab = gsub("10E8V4", "10E8v4", bnab),
                     prediction_performance = est) %>% 
              select(bnab, datatype, prediction_performance), 
            by = c("bnab", "datatype")) %>% 
  mutate(bnab = gsub("_", " + ", bnab),
         bnab = ifelse(bnab == "VRC01-PGDM1400-10E8v4", "VRC01/PGDM1400/10E8v4", bnab))
mc_vars <- all_output %>% 
  group_by(bnab, epsilon, datatype, estimator, simplified) %>% 
  summarize(mc_var = var(est), 
            bias = mean(est - truth),
            pred_perf = mean(prediction_performance), .groups = "drop") %>% 
  pivot_wider(names_from = estimator, values_from = c(mc_var, bias)) %>% 
  mutate(re_nonaug_aug = `mc_var_non-augmented` / mc_var_augmented,
         re_aug_oracle = `mc_var_augmented` / mc_var_oracle,
         re_nonaug_oracle = `mc_var_non-augmented` / mc_var_oracle,
         outcome = ifelse(datatype == "binary", "IC80 < 1", "IC80"),
         percentage = factor(epsilon, levels = c(0.5, 1, 2), labels = c("50%", "100%", "200%")),
         bnab = factor(bnab, levels = c(
           "VRC01", "VRC07-523-LS", "PGT121", "VRC26.25", "PGDM1400",
           "VRC07-523-LS + PGT121", "VRC07-523-LS + VRC26.25", "VRC07-523-LS + PGDM1400", "VRC07-523-LS + 10-1074",
           "VRC07-523-LS + PGT121 + PGDM1400", "VRC01/PGDM1400/10E8v4"
         )))

bias_tib <- all_output %>% 
  group_by(bnab, epsilon, datatype, estimator, simplified) %>% 
  summarize(bias = mean(est - truth), pred_perf = mean(prediction_performance), .groups = "drop") %>% 
  mutate(outcome = ifelse(datatype == "binary", "IC80 < 1", "IC80"),
         percentage = factor(epsilon, levels = c(0.5, 1, 2), labels = c("50%", "100%", "200%")),
         bnab = factor(bnab, levels = c(
           "VRC01", "VRC07-523-LS", "PGT121", "VRC26.25", "PGDM1400",
           "VRC07-523-LS + PGT121", "VRC07-523-LS + VRC26.25", "VRC07-523-LS + PGDM1400", "VRC07-523-LS + 10-1074",
           "VRC07-523-LS + PGT121 + PGDM1400", "VRC01/PGDM1400/10E8v4"
         )))

mc_vars_1 <- mc_vars %>% 
  filter(simplified == 1)
bias_tib_1 <- bias_tib %>% 
  filter(simplified == 1)

mc_vars_0 <- mc_vars %>% 
  filter(simplified == 0)
bias_tib_0 <- bias_tib %>% 
  filter(simplified == 0)


point_size <- 3
bias_ylim <- c(-3e-3, 3e-3)

# ALL PLOTS: SIMPLIFIED DGM
# relative efficiency: non-augmented vs augmented ------------------------------
continuous_rel_eff_plot <- mc_vars_1 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = re_nonaug_aug, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

binary_rel_eff_plot <- mc_vars_1 %>% 
  filter(outcome == "IC80 < 1") %>% 
  ggplot(aes(x = pred_perf, y = re_nonaug_aug, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

lgnd <- get_legend(continuous_rel_eff_plot)
ylab <- get_y_axis(continuous_rel_eff_plot)

full_plot <- plot_grid(
  ggplot() + labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)") + 
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  plot_grid(
    continuous_rel_eff_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
    binary_rel_eff_plot + labs(y = NULL) + theme(legend.position = "none",
                                                           strip.text.x = element_blank()),
    nrow = 2, ncol = 1
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

for (filetype in c("png", "pdf")) {
  ggsave(filename = here::here("R_output", paste0("figure_1.", filetype)),
         plot = full_plot, 
         width = 9.5, height = 5, units = "in")
}

# relative efficiency: non-augmented vs oracle ---------------------------------
continuous_rel_eff_plot_nonaug_oracle <- mc_vars_1 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = re_nonaug_oracle, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring auxiliary sequences vs oracle)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

binary_rel_eff_plot_nonaug_oracle <- mc_vars_1 %>% 
  filter(outcome == "IC80 < 1") %>% 
  ggplot(aes(x = pred_perf, y = re_nonaug_oracle, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring auxiliary sequences vs oracle)", 
       x = "Prediction Performance", color = "bnAb") +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

lgnd <- get_legend(continuous_rel_eff_plot_nonaug_oracle)
ylab <- get_y_axis(continuous_rel_eff_plot_nonaug_oracle)

full_plot_nonaug_oracle <- plot_grid(
  ggplot() + labs(y = "Relative Efficiency (ignoring auxiliary sequences vs oracle)") + 
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  plot_grid(
    continuous_rel_eff_plot_nonaug_oracle + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
    binary_rel_eff_plot_nonaug_oracle + labs(y = NULL) + theme(legend.position = "none",
                                                 strip.text.x = element_blank()),
    nrow = 2, ncol = 1
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

for (filetype in c("png", "pdf")) {
  ggsave(filename = here::here("R_output", paste0("figure_1_nonaug_oracle.", filetype)),
         plot = full_plot_nonaug_oracle, 
         width = 9.5, height = 5, units = "in")
}

# relative efficiency: augmented vs oracle ---------------------------------
continuous_rel_eff_plot_aug_oracle <- mc_vars_1 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = re_aug_oracle, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (predicting auxiliary sequences vs oracle)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

binary_rel_eff_plot_aug_oracle <- mc_vars_1 %>% 
  filter(outcome == "IC80 < 1") %>% 
  ggplot(aes(x = pred_perf, y = re_aug_oracle, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (predicting auxiliary sequences vs oracle)", 
       x = "Prediction Performance", color = "bnAb") +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

lgnd <- get_legend(continuous_rel_eff_plot_aug_oracle)
ylab <- get_y_axis(continuous_rel_eff_plot_aug_oracle)

full_plot_aug_oracle <- plot_grid(
  ggplot() + labs(y = "Relative Efficiency (predicting auxiliary sequences vs oracle)") + 
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  plot_grid(
    continuous_rel_eff_plot_aug_oracle + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
    binary_rel_eff_plot_aug_oracle + labs(y = NULL) + theme(legend.position = "none",
                                                               strip.text.x = element_blank()),
    nrow = 2, ncol = 1
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

for (filetype in c("png", "pdf")) {
  ggsave(filename = here::here("R_output", paste0("figure_1_aug_oracle.", filetype)),
         plot = full_plot_aug_oracle, 
         width = 9.5, height = 5, units = "in")
}

# bias for all estimators ------------------------------------------------------
continuous_bias_plot <- bias_tib_1 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = bias, color = bnab, shape = estimator)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(y = "Estimation bias (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  scale_y_continuous(labels = function(x) sprintf("%g", x), limits = bias_ylim) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

binary_bias_plot <- bias_tib_1 %>% 
  filter(outcome == "IC80 < 1") %>% 
  ggplot(aes(x = pred_perf, y = bias, color = bnab, shape = estimator)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(y = "Estimation bias (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  scale_y_continuous(labels = function(x) sprintf("%g", x), limits = bias_ylim) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

lgnd <- get_legend(continuous_bias_plot)
ylab <- get_y_axis(continuous_bias_plot)

full_bias_plot <- plot_grid(
  ggplot() + labs(y = "Estimation bias (ignoring vs using auxiliary sequences)") + 
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  plot_grid(
    continuous_bias_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
    binary_bias_plot + labs(y = NULL) + theme(legend.position = "none",
                                                 strip.text.x = element_blank()),
    nrow = 2, ncol = 1
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

for (filetype in c("png", "pdf")) {
  ggsave(filename = here::here("R_output", paste0("sim_1_bias.", filetype)),
         plot = full_bias_plot, 
         width = 9.5, height = 5, units = "in")
}

# ALL PLOTS: NON-SIMPLIFIED DGM
# relative efficiency: non-augmented vs augmented ------------------------------
continuous_rel_eff_plot <- mc_vars_0 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = re_nonaug_aug, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

continuous_rel_eff_plot_nonaug_oracle <- mc_vars_0 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = re_nonaug_oracle, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring auxiliary sequences vs oracle)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

continuous_rel_eff_plot_aug_oracle <- mc_vars_0 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = re_aug_oracle, color = bnab)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (predicting auxiliary sequences vs oracle)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

lgnd <- get_legend(continuous_rel_eff_plot)

full_plot <- plot_grid(
  plot_grid(
    continuous_rel_eff_plot_nonaug_oracle + theme(legend.position = "none"),
    continuous_rel_eff_plot_aug_oracle + theme(legend.position = "none"),
    continuous_rel_eff_plot + theme(legend.position = "none"),
    nrow = 1, ncol = 3
  ),
  lgnd, nrow = 1, ncol = 2, rel_widths = c(1, .6)
)

for (filetype in c("png", "pdf")) {
  ggsave(filename = here::here("R_output", paste0("figure_1_unsimplified.", filetype)),
         plot = full_plot, 
         width = 26, height = 6, units = "in")
}

# bias for all estimators ------------------------------------------------------
continuous_bias_plot <- bias_tib_0 %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = bias, color = bnab, shape = estimator)) +
  geom_point(size = point_size) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(y = "Estimation bias (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  scale_y_continuous(labels = function(x) sprintf("%g", x), limits = bias_ylim) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

lgnd <- get_legend(continuous_bias_plot)
ylab <- get_y_axis(continuous_bias_plot)

full_bias_plot <- plot_grid(
  ggplot() + labs(y = "Estimation bias (ignoring vs using auxiliary sequences)") + 
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  continuous_bias_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

for (filetype in c("png", "pdf")) {
  ggsave(filename = here::here("R_output", paste0("sim_1_bias_unsimplified.", filetype)),
         plot = full_bias_plot, 
         width = 9.5, height = 5, units = "in")
}
