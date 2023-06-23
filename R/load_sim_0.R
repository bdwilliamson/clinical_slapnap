library("here")
library("tidyverse")
library("data.table")
library("cowplot")
theme_set(theme_cowplot())

sim_0_results <- readRDS(file = here::here("R_output", "sim_0_results.rds"))

avg_results <- sim_0_results %>% 
  group_by(n) %>% 
  summarize(mdn_mse = median(mse))

point_size <- 3
sim_0_fig <- avg_results %>% 
  ggplot(aes(x = n, y = mdn_mse)) +
  geom_point(size = point_size) +
  labs(x = "n", y = "CV-MSE")

ggsave(filename = here::here("R_output", "sim_0.png"),
       sim_0_fig, width = 5, height = 5, units = "in")
