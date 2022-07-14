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
         bnab = ifelse(bnab == "VRC01-PGDM1400-10E8v4", "VRC01/PGDM1400/10E8v4", bnab))
mc_vars <- all_output %>% 
  group_by(bnab, epsilon, datatype, augmented) %>% 
  summarize(mc_var = var(est), 
            pred_perf = mean(prediction_performance), .groups = "drop") %>% 
  pivot_wider(names_from = augmented, values_from = mc_var,
              names_prefix = "aug_") %>% 
  mutate(relative_efficiency = aug_FALSE / aug_TRUE,
         outcome = ifelse(datatype == "binary", "IC80 < 1", "IC80"),
         percentage = factor(epsilon, levels = c(0.5, 1, 2), labels = c("50%", "100%", "200%")),
         bnab = factor(bnab, levels = c(
           "VRC01", "VRC07-523-LS", "PGT121", "VRC26.25", "PGDM1400",
           "VRC07-523-LS + PGT121", "VRC07-523-LS + VRC26.25", "VRC07-523-LS + PGDM1400", "VRC07-523-LS + 10-1074",
           "VRC07-523-LS + PGT121 + PGDM1400", "VRC01/PGDM1400/10E8v4"
         )))

continuous_rel_eff_plot <- mc_vars %>% 
  filter(outcome == "IC80") %>% 
  ggplot(aes(x = pred_perf, y = relative_efficiency, color = bnab)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)", 
       x = "Prediction Performance", color = "bnAb") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
                                         breaks = NULL, labels = NULL)) +
  facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both)

binary_rel_eff_plot <- mc_vars %>% 
  filter(outcome == "IC80 < 1") %>% 
  ggplot(aes(x = pred_perf, y = relative_efficiency, color = bnab)) +
  geom_point() +
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
    binary_rel_eff_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none",
                                                           strip.text.x = element_blank()),
    nrow = 2, ncol = 1
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

# rel_eff_plot <- mc_vars %>%
#   ggplot(aes(x = pred_perf, y = relative_efficiency, color = bnab)) +
#   geom_point() +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
#   labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)",
#        x = "Prediction Performance", color = "bnAb") +
#   scale_x_continuous(sec.axis = sec_axis(~ ., name = "Percentage Increase in Available Viruses",
#                      breaks = NULL, labels = NULL)) +
#   facet_grid(rows = vars(outcome), cols = vars(percentage), labeller = label_both,
#              scales = "free")

ggsave(filename = here::here("R_output", "sim_1_rel_eff.png"),
       plot = full_plot, 
       width = 9, height = 5, units = "in")
