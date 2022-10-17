#!/usr/local/bin/Rscript
library("here")
library("tidyverse")
library("data.table")
library("cowplot")
theme_set(theme_cowplot())

# read in results of the simulation, create a table
all_output_files <- list.files(here::here("R_output", "simulation_2"))
all_output <- tibble::tibble(data.table::rbindlist(
  lapply(as.list(here::here("R_output", "simulation_2", all_output_files)), readRDS)
))
pe_0_labels <- sort(round(unique(all_output$pe_0), 3))
results <- all_output %>% 
  group_by(analysis, n, position, pe_0) %>% 
  summarize(power = mean(reject), .groups = "drop") %>% 
  mutate(analysis = factor(analysis), n = factor(n), 
         `PE(S230 = 0)` = factor(round(pe_0, 3),
                                 levels = pe_0_labels))

power_plot <- results %>% 
  arrange(`PE(S230 = 0)`) %>% 
  ggplot(aes(x = n, y = power, color = analysis, group = analysis, shape = analysis)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(y = "Empirical power", x = "n", color = "Analysis", shape = "Analysis") +
  facet_grid(cols = vars(`PE(S230 = 0)`), labeller = label_both) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_color_viridis_d(begin = 0, end = 0.5)

ggsave(filename = here::here("R_output", paste0("sim_2.png")),
       plot = power_plot,
       width = 9, height = 5, units = "in")
ggsave(filename = here::here("R_output", paste0("figure_3.pdf")),
       plot = power_plot,
       width = 9, height = 5, units = "in")
