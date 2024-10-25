library(tidyverse)
library(gridExtra)
library(cowplot)
source("4_plot_helper_functions.R")

dwbin <- dwbin <- .6
axis_title_size <- 12
axis_text_size <- 11
strip_text_size <- 12
subtitle_size <- 12
legend_size <- 13
legend_icon_size <- 2
legend_dot_size <- 3.5

results_binary <- readRDS("results/results_irgtt_binary.rds") %>%
  filter(treatment_arm_clusters > 12) %>% 
  mutate(style = "irgtt")
  
  
bin_complex <- summarize_logit(results_binary,
                       version = "complex_covars")
plot1 <- make_grid_plot(bin_complex, continuous = F)
ggsave(plot = plot1, dpi = 300,
       width = 10, height = 10, 
       filename = "irgtt_complex_binary.pdf")


bin_main_terms <- summarize_logit(results_binary,
                               version = "main_terms")
ggsave(plot = make_grid_plot(bin_main_terms, continuous = F), dpi = 300,
       width = 10, height = 10, 
       filename = "irgtt_main_terms_binary.pdf")


bin_overfit <- summarize_logit(results_binary,
                                  version = "overfit")
ggsave(plot = make_grid_plot(bin_overfit, continuous = F), dpi = 300,
       width = 10, height = 10, 
       filename = "irgtt_overfit_binary.pdf")



################################
strip_text_size <- 9
results_binary_ub <- readRDS("results/results_irgtt_binary_unbalanced.rds") %>%
  mutate(style = "unbalanced")

bin_complex_ub <- summarize_logit(results_binary_ub,
                               version = "complex_covars")
plot1_ub <- make_grid_plot(bin_complex_ub, continuous = F)
ggsave(plot = plot1_ub, dpi = 300,
       width = 10, height = 10, 
       filename = "irgtt_complex_binary_ub.pdf")


bin_main_terms_ub <- summarize_logit(results_binary_ub,
                                  version = "main_terms")
ggsave(plot = make_grid_plot(bin_main_terms_ub, continuous = F), dpi = 300,
       width = 10, height = 10, 
       filename = "irgtt_main_terms_binary_ub.pdf")


bin_overfit_ub <- summarize_logit(results_binary_ub,
                               version = "overfit")
ggsave(plot = make_grid_plot(bin_overfit_ub, continuous = F), dpi = 300,
       width = 10, height = 10, 
       filename = "irgtt_overfit_binary_ub.pdf")






