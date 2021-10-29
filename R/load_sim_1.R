#!/usr/local/bin/Rscript
library("here")
library("tidyverse")
library("data.table")
library("cowplot")
theme_set(theme_cowplot())

# read in results of the simulation, create a table
all_output_files <- list.files(here::here("R_output", "simulation_1"))
all_output <- tibble::tibble(data.table::rbindlist(
  lapply(as.list(here::here("R_output", "simulation_1", all_output_files)), readRDS)
))
mc_vars <- all_output %>% 
  group_by(bnab, epsilon, datatype, augmented) %>% 
  summarize(mc_var = var(est), .groups = "drop") %>% 
  pivot_wider(names_from = augmented, values_from = mc_var,
              names_prefix = "aug_") %>% 
  mutate(relative_efficiency = aug_FALSE / aug_TRUE,
         outcome = ifelse(datatype == "binary", "IC80 < 1", "IC80"))


rel_eff_plot <- mc_vars %>% 
  ggplot(aes(x = factor(epsilon), y = relative_efficiency, color = bnab)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylim(c(0.75, 2.5)) +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)", x = "Fraction of auxiliary sequences", color = "bnAb") +
  facet_grid(rows = vars(outcome), labeller = label_both)

ggsave(filename = here::here("R_output", "sim_1_rel_eff.png"),
       plot = rel_eff_plot, 
       width = 20, height = 15, units = "cm")
