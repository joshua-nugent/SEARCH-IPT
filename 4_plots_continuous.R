library(tidyverse)
library(gridExtra)
library(cowplot)
source("4_plot_helper_functions.R")

dwbin <- dwcont <- .6
axis_title_size <- 8
axis_text_size <- 7
strip_text_size <- 8
subtitle_size <- 8
legend_size <- 9
legend_icon_size <- 1
legend_dot_size <- 2

###################################
results_continuous <- readRDS("results/results_irgtt.rds") %>%
  filter(treatment_arm_clusters > 12) %>% 
  mutate(style = "irgtt")

cont_complex <- summarize_continuous(results_continuous, version = "complex_covars")
cont_complex %>% filter(treatment_arm_cluster_size == 5, treatment_arm_clusters == 50, ranef_sd == .23,t_eff >0) %>% 
  select(contains("power"))
mean(cont_complex$tmle_cov[cont_complex$t_eff > 0])
mean(cont_complex$tmle_TIE[cont_complex$t_eff==0])
ggsave(device = "pdf", plot = make_grid_plot(cont_complex),
       width = 10, height = 10, dpi = 300,
       filename = "irgtt_complex_continuous.pdf")


cont_main_terms <- summarize_continuous(results_continuous, version = "main_terms")
ggsave(device = "pdf", plot = make_grid_plot(cont_main_terms),
       width = 10, height = 10, dpi = 300,
       filename = "irgtt_main_terms_continuous.pdf")

cont_overfit <- summarize_continuous(results_continuous, version = "overfit")
ggsave(device = "pdf", plot = make_grid_plot(cont_overfit),       
       width = 10, height = 10, dpi = 300,
       filename = "irgtt_overfit_continuous.pdf")


###################################
results_continuous_ub <- readRDS("results/results_irgtt_unbalanced.rds") %>%
  mutate(style = "unbalanced")

cont_complex_ub <- summarize_continuous(results_continuous_ub, version = "complex_covars")
p_complex_cont_ub <- make_grid_plot(cont_complex_ub)
ggsave(device = "pdf", plot = p_complex_cont_ub,
       width = 10, height = 10, dpi = 300,
       filename = "irgtt_complex_continuous_ub.pdf")

cont_main_terms_ub <- summarize_continuous(results_continuous_ub, version = "main_terms")
ggsave(device = "pdf", plot = make_grid_plot(cont_main_terms_ub),
       width = 10, height = 10, dpi = 300,
       filename = "irgtt_main_terms_continuous_ub.pdf")

cont_overfit_ub <- summarize_continuous(results_continuous_ub, version = "overfit")
ggsave(device = "pdf", plot = make_grid_plot(cont_overfit_ub),
       width = 10, height = 10, dpi = 300,
       filename = "irgtt_overfit_continuous_ub.pdf")























