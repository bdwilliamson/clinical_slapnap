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
         sq_width = ifelse(sq_width == 0, NA, sq_width),
         est_ic80 = factor(ifelse(outcome == "ic80", as.numeric(est < 0), as.numeric(est > 0.5))))
ci_widths <- all_output %>%
  select(-est, -starts_with("ci"), -width) %>%
  pivot_wider(names_from = augmented, values_from = c(sq_width, est_ic80),
              names_prefix = "aug_") %>%
  mutate(relative_efficiency = sq_width_aug_FALSE / sq_width_aug_TRUE,
         est_ic80 = est_ic80_aug_FALSE,
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
           outcome == "IC80 < 1" & prediction_performance <= 0.65 ~ "1",
           outcome == "IC80 < 1" & prediction_performance > 0.65 ~ "2",
         ), levels = c("1", "2"), labels = c("[0.5, 0.65]", "(0.65, 1]")),
         `R-squared` = factor(case_when(
           outcome == "IC80" & prediction_performance <= 0.32 ~ "1",
           outcome == "IC80" & prediction_performance > 0.32 ~ "2",
         ), levels = c("1", "2"), labels = c("[0, 0.32]", "(0.32, 1]")),
         n_lanl = n_overall - n_catnap,
         excess_prop_catnap = (n_lanl - n_catnap) / n_catnap * 100,
         excess_prop_lanl = (n_lanl - n_catnap) / n_lanl,
         efficiency_bound = n_overall / n_catnap,
         bounded_relative_efficiency = pmin(relative_efficiency, efficiency_bound))

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
  ggplot(aes(x = excess_prop_lanl, y = relative_efficiency, color = bnab, shape = country)) +
  geom_point(size = 2.5) +
  scale_shape_manual(values = c(49:56)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)",
       x = "Proportion of Additional Information in LANL", color = "bnAb",
       shape = "Country") +
  facet_grid(rows = vars(outcome), labeller = label_both)

binary_rel_eff_plot <- ci_widths %>%
  filter(outcome == "IC80 < 1") %>%
  ggplot(aes(x = excess_prop_lanl, y = relative_efficiency, color = bnab, shape = country)) +
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


ci_widths %>%
  filter(outcome == "IC80", country == "United States") %>% 
  select(bnab, n_catnap, n_overall, prediction_performance, relative_efficiency)

all_output %>% 
  filter(outcome == "ic80", country == "US") %>% 
  select(bnab, n_catnap, n_overall, ci_ll, ci_ul, sq_width) %>% 
  print(n = Inf)

all_output %>% 
  filter(outcome == "ic80", country == "US") %>% 
  select(bnab, n_catnap, n_overall, prediction_performance, augmented, est, sq_width) %>% 
  mutate(transformed_est = 10 ^ est) %>% 
  print(n = Inf)

ci_widths %>%
  filter(outcome == "IC80", country == "China") %>% 
  select(bnab, n_catnap, n_overall, prediction_performance, relative_efficiency, bounded_relative_efficiency)

ci_widths %>% 
  filter(n_catnap > 50) %>% 
  select(bnab, country, n_catnap, n_overall, prediction_performance, relative_efficiency, bounded_relative_efficiency) %>% 
  print(n = Inf)
  
# based on the bounds
dodge_width <- 0.01
dodge_height <- 0.75
continuous_rel_eff_plot_bounded <- ci_widths %>%
  filter(outcome == "IC80") %>%
  ggplot(aes(x = excess_prop_lanl, y = bounded_relative_efficiency, color = bnab, shape = country)) +
  geom_point(size = 2.5, position = position_jitter(height = dodge_height)) +
  scale_shape_manual(values = c(49:56)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)",
       x = "Proportion of Additional Information in LANL", color = "bnAb",
       shape = "Country") +
  facet_grid(rows = vars(outcome), labeller = label_both)

binary_rel_eff_plot_bounded <- ci_widths %>%
  filter(outcome == "IC80 < 1") %>%
  ggplot(aes(x = excess_prop_lanl, y = bounded_relative_efficiency, color = bnab, shape = country)) +
  geom_point(size = 2.5, position = position_jitter(width = dodge_width, height = dodge_height)) +
  scale_shape_manual(values = c(49:56)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)",
       x = "Proportion of Additional Information in LANL", color = "bnAb",
       shape = "Country") +
  facet_grid(rows = vars(outcome), labeller = label_both)

lgnd_b <- get_legend(continuous_rel_eff_plot_bounded)
ylab_b <- get_y_axis(continuous_rel_eff_plot_bounded)

full_plot_b <- plot_grid(
  ggplot() + labs(y = "Relative Efficiency (ignoring vs using auxiliary sequences)") +
    theme(axis.line.x = element_blank(), axis.line.y = element_blank()),
  plot_grid(
    continuous_rel_eff_plot_bounded + labs(x = NULL, y = NULL) + theme(legend.position = "none"),
    binary_rel_eff_plot_bounded + labs( y = NULL) + theme(legend.position = "none"),
    nrow = 2, ncol = 1, rel_heights = c(0.9, 1)
  ),
  lgnd, nrow = 1, ncol = 3, rel_widths = c(.05, 1, .6)
)

ggsave(filename = here::here("R_output", "sim_1b_rel_eff_bounded.png"),
       plot = full_plot_b,
       width =  11.5, height = 5, units = "in")
