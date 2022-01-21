#!/usr/local/bin/Rscript
library("here")
library("tidyverse")
library("data.table")
library("cowplot")
theme_set(theme_cowplot())

summary_statistics <- readRDS(here::here("R_output", "summary_statistics_simulation_1.rds"))

# read in results of the simulation, create a table
all_output_files <- list.files(here::here("R_output", "simulation_1b"))
all_output <- tibble::tibble(data.table::rbindlist(
  lapply(as.list(here::here("R_output", "simulation_1b", all_output_files)), readRDS)
)) %>% 
  left_join(summary_statistics %>% 
              mutate(datatype = ifelse(outcome == "ic80", "continuous", "binary"),
                     bnab = toupper(bnab),
                     bnab = gsub("10E8V4", "10E8v4", bnab),
                     prediction_performance = est) %>% 
              select(bnab, outcome, datatype, prediction_performance), 
            by = c("bnab", "outcome")) %>% 
  mutate(bnab = gsub("_", " + ", bnab),
         bnab = ifelse(bnab == "VRC01-PGDM1400-10E8v4", "VRC01/PGDM1400/10E8v4", bnab),
         epsilon = n_overall / n_catnap)

# Bin epsilon into reasonable groups, like in the simulated data
quantile(all_output$epsilon, seq(0, 1, .33))

# Bin and average within the bins
results <- all_output %>% 
  mutate(percentage = factor(
    case_when(
      epsilon < 8 ~ 1,
      8 <= epsilon & epsilon <= 11 ~ 2,
      11 < epsilon ~ 3
    ), labels = c("< 800%", "800--1100%", "> 1100%"))) %>% 
  select(-est, -ci_ll, -ci_ul) %>% 
  pivot_wider(names_from = augmented, values_from = width,
              names_prefix = "aug_") %>% 
  mutate(relative_efficiency = aug_FALSE / aug_TRUE) %>% 
  group_by(bnab, outcome, prediction_performance, percentage) %>% 
  summarize(relative_efficiency = mean(relative_efficiency), .groups = "drop")

# plot: ci width (y-axis) vs n(catnap) / n(overall), colored by bnab, shape = augmented, panelled by outcome?
continuous_rel_eff_plot <- results %>% 
  filter(outcome == "ic80") %>% 
  ggplot(aes(x = prediction_performance, y = relative_efficiency, color = bnab)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), 
             labeller = label_bquote(rows = Outcome:~IC[80]))

binary_rel_eff_plot <- results %>% 
  filter(outcome == "sens") %>% 
  ggplot(aes(x = prediction_performance, y = relative_efficiency, color = bnab)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  facet_grid(rows = vars(outcome), cols = vars(percentage), 
             labeller = label_bquote(rows = Outcome:~IC[80] < 1))

lgnd <- get_legend(continuous_rel_eff_plot)
ylab <- get_y_axis(continuous_rel_eff_plot)

full_plot <- plot_grid(
  ggplot() + labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)") + 
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  plot_grid(
    continuous_rel_eff_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
    binary_rel_eff_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none",
                                                           strip.text.x = element_blank()),
    nrow = 2, ncol = 1
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

# figure for the main paper
ggsave(filename = here::here("R_output", "sim_1b_rel_eff.png"),
       plot = full_plot, 
       width = 9, height = 5, units = "in")

# table for the supplement
all_output %>% 
  select(-datatype) %>% 
  write_csv(path = here::here("R_output", "sim_1b_results.csv"))