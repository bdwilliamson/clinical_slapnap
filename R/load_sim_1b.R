#!/usr/local/bin/Rscript
library("here")
library("tidyverse")
library("data.table")
library("cowplot")
theme_set(theme_cowplot())

summary_statistics <- readRDS(here::here("R_output", "summary_statistics_simulation_1.rds"))
summary_stats <- summary_statistics %>%
  mutate(datatype = ifelse(outcome == "ic80", "continuous", "binary"),
         bnab = toupper(bnab),
         bnab = gsub("10E8V4", "10E8v4", bnab),
         prediction_performance = est) %>%
  select(bnab, outcome, datatype, prediction_performance)
# read in results of the simulation, create a table
all_output_files <- list.files(here::here("R_output", "simulation_1b"))
all_output <- tibble::tibble(data.table::rbindlist(
  lapply(as.list(here::here("R_output", "simulation_1b", all_output_files)), readRDS)
)) %>%
  left_join(summary_stats, by = c("bnab", "outcome")) %>%
  mutate(bnab = gsub("_", " + ", bnab),
         bnab = ifelse(bnab == "VRC01-PGDM1400-10E8v4", "VRC01/PGDM1400/10E8v4", bnab),
         width = ifelse(width == 0, NA, width),
         sq_width = ifelse(sq_width == 0, NA, sq_width))
ci_widths <- all_output %>%
  select(-est, -starts_with("ci"), -width) %>%
  # group_by(bnab, country, n_catnap, n_overall, datatype, augmented) %>%
  pivot_wider(names_from = augmented, values_from = sq_width,
              names_prefix = "aug_") %>%
  mutate(relative_efficiency = aug_FALSE / aug_TRUE,
         outcome = ifelse(datatype == "binary", "IC80 < 1", "IC80"),
         # percentage = factor(epsilon, levels = c(0.5, 1, 2), labels = c("50%", "100%", "200%")),
         bnab = factor(bnab, levels = c(
           "VRC01", "VRC07-523-LS", "PGT121", "VRC26.25", "PGDM1400",
           "VRC07-523-LS + PGT121", "VRC07-523-LS + VRC26.25", "VRC07-523-LS + PGDM1400", "VRC07-523-LS + 10-1074",
           "VRC07-523-LS + PGT121 + PGDM1400", "VRC01/PGDM1400/10E8v4"
         )),
         country = factor(country, levels = c(
           "CN", "DE", "KE", "MW", "TZ", "UG", "US", "ZA"
         ), labels = c(
           "China", "Germany", "Kenya", "Malawi", "Tanzania", "Uganda", "United States", "South Africa"
         )),
         `AUC` = factor(case_when(
           outcome == "IC80 < 1" & prediction_performance <= 0.7 ~ "1",
           outcome == "IC80 < 1" & prediction_performance > 0.7 & prediction_performance <= 0.8 ~ "2",
           outcome == "IC80 < 1" & prediction_performance > 0.8 ~ "3",
         ), levels = c("1", "2", "3"), labels = c("[0.5, 0.7]", "(0.7, 0.8]", "(0.8, 1]")),
         `R-squared` = factor(case_when(
           outcome == "IC80" & prediction_performance <= 0.2 ~ "1",
           outcome == "IC80" & prediction_performance > 0.2 & prediction_performance <= 0.4 ~ "2",
           outcome == "IC80" & prediction_performance > 0.4 ~ "3",
         ), levels = c("1", "2", "3"), labels = c("[0, 0.2]", "(0.2, 0.4]", "[0.4, 1]")),
         excess_prop = (n_overall - n_catnap) / n_catnap * 100)

# Distribution of numbers in CATNAP, LANL
numbers_only <- all_output %>% 
  select(bnab, country, n_catnap, n_overall) %>% 
  group_by(bnab, country, n_catnap, n_overall) %>% 
  slice(1) %>% 
  pivot_longer(cols = starts_with("n_"), names_to = "source", values_to = "n") %>% 
  mutate(source = factor(gsub("n_", "", source), levels = c("overall", "catnap"), labels = c("LANL", "CATNAP")),
         bnab_country = factor(paste0(bnab, "_", country)))
dist_catnap <- numbers_only %>% 
  filter(source == "CATNAP") %>% 
  ggplot(aes(x = country, y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "Country", y = "Number of sequences", fill = "Data source") +
  facet_wrap(~ bnab)

# Relative efficiency
continuous_rel_eff_plot <- ci_widths %>%
  filter(outcome == "IC80") %>%
  ggplot(aes(x = excess_prop, y = relative_efficiency, color = bnab, shape = country)) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(49:56)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)",
       x = "Proportion of Additional Information in LANL", color = "bnAb",
       shape = "Country") +
  facet_grid(rows = vars(outcome), labeller = label_both)

binary_rel_eff_plot <- ci_widths %>%
  filter(outcome == "IC80 < 1") %>%
  ggplot(aes(x = excess_prop, y = relative_efficiency, color = bnab, shape = country)) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(49:56)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)",
       x = "Proportion of Additional Information in LANL", color = "bnAb",
       shape = "Country") +
  facet_grid(rows = vars(outcome), labeller = label_both)

lgnd <- get_legend(continuous_rel_eff_plot)
ylab <- get_y_axis(continuous_rel_eff_plot)

full_plot <- plot_grid(
  ggplot() + labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)") +
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  plot_grid(
    continuous_rel_eff_plot + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
    binary_rel_eff_plot + labs( y = NULL) + theme(legend.position = "none"),
    nrow = 2, ncol = 1, rel_heights = c(0.9, 1)
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

ggsave(filename = here::here("R_output", "sim_1b_rel_eff.png"),
       plot = full_plot,
       width =  11.5, height = 5, units = "in")
