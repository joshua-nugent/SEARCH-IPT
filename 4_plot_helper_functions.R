library(tidyverse)

make_grid_plot <- function(plotdat, continuous = T){
  if(continuous){
    grid.arrange(generate_bias_plot(plotdat),
                 generate_power_plot(plotdat),
                 generate_MSE_plot(plotdat),
                 generate_TIE_plot(plotdat),
                 generate_coverage_plot(plotdat),
                 generate_legend_plot(plotdat),
                 ncol = 2)
  } else {
    grid.arrange(generate_power_plot(plotdat),
                 generate_TIE_plot(plotdat),
                 generate_var_ratio_plot(plotdat),
                 generate_legend_plot(plotdat),
                 ncol = 2) 
  }
}


summarize_logit <- function(results, version = "complex_covars"){
  results %>% filter(link_function == "logit",
                     sim_type == version,
                     !is.na(glmm_se),
                     !is.na(glmm_pe),
                     !is.infinite(glmm_pe),
                     !is.infinite(gee_se),
                     !is.nan(glmm_lo),
                     !is.nan(glmm_hi)
  ) %>%
    mutate(ranef_sd = round(ranef_sd, digits = 2)) %>% 
    mutate(across(contains("ranef_sd"), as.factor)) %>% 
    mutate(across(contains("t_eff"), as.factor)) %>% 
    group_by(style,
             link_function,
             control_arm_irgtt,
             control_arm_no_ranef,
             control_arm_clusters,
             treatment_arm_clusters,
             control_arm_cluster_size,
             treatment_arm_cluster_size,
             sizes,
             t_eff,
             ranef_sd,
             sim_type, nsims
    ) %>%
    summarize(n = n(),
              
              tmle_mean_se = mean(ltmle_se, na.rm = T),
              glmm_mean_se = mean(glmm_se, na.rm = T),
              gee_mean_se = mean(gee_se, na.rm = T),

              tmle_pe_sd = sd(log(ltmle_pe), na.rm = T),
              glmm_pe_sd = sd(log(glmm_pe), na.rm = T),
              gee_pe_sd = sd(log(gee_pe), na.rm = T),

              tmle_var_ratio = tmle_mean_se / tmle_pe_sd,
              glmm_var_ratio = glmm_mean_se / glmm_pe_sd,
              gee_var_ratio = gee_mean_se / gee_pe_sd,

              tmle_TIE = mean(ltmle_lo > 1 | ltmle_hi < 1, na.rm = T),
              glmm_TIE = mean(glmm_lo > 1 | glmm_hi < 1, na.rm = T),
              gee_TIE = mean(gee_lo > 1 | gee_hi < 1, na.rm = T),

              tmle_power = mean(ltmle_lo > 1 | ltmle_hi < 1, na.rm = T),
              glmm_power = mean(glmm_lo > 1 | glmm_hi < 1, na.rm = T),
              gee_power = mean(gee_lo > 1 | gee_hi < 1, na.rm = T),

              tmle_cise_hi = tmle_mean_se + 1.96*sqrt(var(ltmle_se) / (4*n*tmle_mean_se^2)),
              glmm_cise_hi = glmm_mean_se + 1.96*sqrt(var(glmm_se) / (4*n*glmm_mean_se^2)),
              gee_cise_hi = gee_mean_se + 1.96*sqrt(var(gee_se) / (4*n*gee_mean_se^2)),

              tmle_cise_lo = tmle_mean_se - 1.96*sqrt(var(ltmle_se) / (4*n*tmle_mean_se^2)),
              glmm_cise_lo = glmm_mean_se - 1.96*sqrt(var(glmm_se) / (4*n*glmm_mean_se^2)),
              gee_cise_lo = gee_mean_se - 1.96*sqrt(var(gee_se) / (4*n*gee_mean_se^2)),
    ) %>% 
    mutate(`Control arm modeled as` = ifelse(control_arm_irgtt == 1, "unclustered", "clustered"),
           `Treatment effect` = t_eff) %>%
    ungroup()
}


summarize_continuous <- function(results, version = "complex_covars"){
  output <- results %>% filter(link_function == "identity",
                               sim_type == version,
                               !is.na(glmm_se),
                               !is.na(glmm_pe),
                               !is.infinite(glmm_pe),
                               !is.infinite(gee_se),
                               !is.infinite(gee_hi),
                               !is.nan(glmm_lo),
                               !is.nan(glmm_hi),
                               !is.nan(gee_lo),
                               !is.nan(gee_hi),
                               glmm_se < 5,
                               gee_se < 5,
                               abs(glmm_pe) < 10) %>%
    mutate(ranef_sd = round(ranef_sd, digits = 2)) %>% 
    mutate(across(contains("ranef_sd"), as.factor)) %>% 
    #mutate(across(contains("t_eff"), as.factor)) %>% 
    group_by(style,
             link_function,
             control_arm_irgtt,
             control_arm_no_ranef,
             control_arm_clusters,
             treatment_arm_clusters,
             control_arm_cluster_size,
             treatment_arm_cluster_size,
             sizes,
             t_eff,
             ranef_sd,
             #v, #true_marg_mean_diff,
             sim_type,
             nsims
    ) %>%
    summarize(n = n(),
              
              tmle_mean_se = mean(ltmle_se, na.rm = T),
              glmm_mean_se = mean(glmm_se, na.rm = T),
              gee_mean_se = mean(gee_se, na.rm = T),

              tmle_pe_sd = sd(ltmle_pe, na.rm = T),
              glmm_pe_sd = sd(glmm_pe, na.rm = T),
              gee_pe_sd = sd(gee_pe, na.rm = T),

              tmle_var_ratio = tmle_mean_se / tmle_pe_sd,
              glmm_var_ratio = glmm_mean_se / glmm_pe_sd,
              gee_var_ratio = gee_mean_se / gee_pe_sd,

              tmle_mse = mean((ltmle_pe - t_eff)^2, na.rm = T),
              glmm_mse = mean((glmm_pe - t_eff)^2, na.rm = T),
              gee_mse = mean((gee_pe - t_eff)^2, na.rm = T),

              tmle_semse = sd((ltmle_pe - t_eff)^2, na.rm = T) / sqrt(n),
              glmm_semse = sd((glmm_pe - t_eff)^2, na.rm = T) / sqrt(n),
              gee_semse = sd((gee_pe - t_eff)^2, na.rm = T) / sqrt(n),

              tmle_himse = tmle_mse + 1.96*tmle_semse,
              glmm_himse = glmm_mse + 1.96*glmm_semse,
              gee_himse = gee_mse + 1.96*gee_semse,

              tmle_lomse = tmle_mse - 1.96*tmle_semse,
              glmm_lomse = glmm_mse - 1.96*glmm_semse,
              gee_lomse = gee_mse - 1.96*gee_semse,

              tmle_TIE = mean(ltmle_lo > 0 | ltmle_hi < 0),
              glmm_TIE = mean(glmm_lo > 0 | glmm_hi < 0),
              gee_TIE = mean(gee_lo > 0 | gee_hi < 0),

              tmle_power = mean(ltmle_lo > 0 | ltmle_hi < 0, na.rm = T),
              glmm_power = mean(glmm_lo > 0 | glmm_hi < 0, na.rm = T),
              gee_power = mean(gee_lo > 0 | gee_hi < 0, na.rm = T),

              tmle_cov = mean(ltmle_lo < t_eff & ltmle_hi > t_eff, na.rm = T),
              glmm_cov = mean(glmm_lo < t_eff & glmm_hi > t_eff, na.rm = T),
              gee_cov = mean(gee_lo < t_eff & gee_hi > t_eff, na.rm = T),

              tmle_bias = mean(ltmle_pe - t_eff, na.rm = T),
              glmm_bias = mean(glmm_pe - t_eff, na.rm = T),
              gee_bias = mean(gee_pe - t_eff, na.rm = T),

              tmle_cibias_hi = tmle_bias + 1.96*tmle_pe_sd/sqrt(n),
              glmm_cibias_hi = glmm_bias + 1.96*glmm_pe_sd/sqrt(n),
              gee_cibias_hi = gee_bias + 1.96*gee_pe_sd/sqrt(n),

              tmle_cibias_lo = tmle_bias - 1.96*tmle_pe_sd/sqrt(n),
              glmm_cibias_lo = glmm_bias - 1.96*glmm_pe_sd/sqrt(n),
              gee_cibias_lo = gee_bias - 1.96*gee_pe_sd/sqrt(n),

              # variance of SE estimates from the models
              tmle_cise_hi = tmle_mean_se + 1.96*sqrt(var(ltmle_se) / (4*n*tmle_mean_se^2)),
              glmm_cise_hi = glmm_mean_se + 1.96*sqrt(var(glmm_se) / (4*n*glmm_mean_se^2)),
              gee_cise_hi = gee_mean_se + 1.96*sqrt(var(gee_se) / (4*n*gee_mean_se^2)),

              tmle_cise_lo = tmle_mean_se - 1.96*sqrt(var(ltmle_se) / (4*n*tmle_mean_se^2)),
              glmm_cise_lo = glmm_mean_se - 1.96*sqrt(var(glmm_se) / (4*n*glmm_mean_se^2)),
              gee_cise_lo = gee_mean_se - 1.96*sqrt(var(gee_se) / (4*n*gee_mean_se^2)),
    ) %>% mutate(across(contains("t_eff"), as.factor)) %>% 
    mutate(#`Control arm modeled as` = ifelse(control_arm_irgtt == 1, "unclustered", "clustered"),
      `Treatment effect` = t_eff) %>%
    mutate(across(contains("Treatment effect"), as.factor)) %>% 
    ungroup()
  return(output)
}

#results$sizes[[1]]
#cont <- summarize_continuous(results)

get_facet <- function(dat, style, type = "bias"){
  
  if(type %in% c("bias", "TIE", "coverage", "ratio", "power")){
    scale_setting <- NULL
    if(type == "power" & dat$link_function[1] == "logit"){
      scale_setting <- "free_y"
    }
  } else {
    scale_setting <- "free"
  }
  
  
  if(style == "balanced"){
    dat <- dat %>% rename(cluster_size = control_arm_cluster_size,
                          clusters_per_arm = control_arm_clusters)
    facet <- facet_grid(clusters_per_arm ~ cluster_size, scales = scale_setting,
                        labeller = "label_both")
  } else if(style == "irgtt"){
    dat <- dat %>% rename(`N<sub>k</sub>` = control_arm_cluster_size,
                          K = control_arm_clusters)
    facet <- facet_grid(K ~ `N<sub>k</sub>`, scales = scale_setting,
                        labeller = "label_both")
    # dat <- dat %>% rename(cluster_size = control_arm_cluster_size,
    #                       clusters_per_arm = control_arm_clusters)
    # facet <- facet_grid(clusters_per_arm ~ cluster_size, scales = scale_setting,
    #                     labeller = "label_both")
  } else if(style == "unbalanced"){
    versions <- paste0(dat$control_arm_clusters * dat$control_arm_cluster_size,
                       " control arm individuals,<br>", 
                       dat$treatment_arm_clusters," treatment arm clusters of size ", 
                       dat$treatment_arm_cluster_size)
    dat <- cbind.data.frame(dat, versions = versions)
    facet <- facet_wrap(~versions, scales = scale_setting)
  } else if(style == "irgtt_unbalanced"){
    versions <- paste0("control clusters = ", dat$control_arm_clusters,
                       ", control cluster size = ", dat$control_arm_cluster_size,
                       ",<br>treatment clusters = ", dat$treatment_arm_clusters,
                       ", treatment cluster size = ", dat$treatment_arm_cluster_size)
    dat <- cbind.data.frame(dat, versions = versions)
    facet <- facet_wrap(~as.character(versions), scales = scale_setting)
  }
  return(list(dat = dat, facet = facet))
}

##########################################
# bias
##########################################
generate_bias_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "bias")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("_bias"), names_to = "method", values_to = "bias") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  plotdat_hi <- dat %>% #filter(t_eff != 0) %>%
    pivot_longer(cols = contains("_cibias_hi"), names_to = "method", values_to = "cibias_hi") %>%
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  plotdat_lo <- dat %>% #filter(t_eff != 0) %>%
    pivot_longer(cols = contains("_cibias_lo"), names_to = "method", values_to = "cibias_lo") %>%
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  pdci <- inner_join(plotdat_lo, plotdat_hi)
  
  if(exclude_clustered){
    pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot(data = plotdat,
         aes(x = ranef_sd,
             group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`)),
         color = method) +
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = .5) +
    geom_point(aes(y = bias,
                   color = method,
                   shape = `Treatment effect`),
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(data = pdci,
                  aes(x = ranef_sd,
                      color = method,
                      group = interaction(`Treatment effect`, method, ranef_sd, `Control arm modeled as`),
                      ymin = cibias_lo, ymax = cibias_hi,
                      linetype = `Control arm modeled as`),
                  width = .3, position = position_dodge(width = dodge_width)) +
    facet +
    theme(legend.position = "none",
          title = element_text(size = subtitle_size),
          axis.text = element_text(size = axis_text_size),
          #strip.text = element_text(size = strip_text_size),
          strip.text = ggtext::element_markdown(size = strip_text_size),
          axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept", title = "Bias", y = NULL)
}

#results <- readRDS("results/results_balanced.rds") %>% mutate(style = "balanced")
#results <- readRDS("results/results_irgtt.rds") %>% mutate(style = "irgtt")
#results <- readRDS("results/results_unbalanced.rds") %>% mutate(style = "unbalanced")
#results <- readRDS("results/results_IRGTT_unbalanced.rds") %>% mutate(style = "irgtt_unbalanced")
#cont <- summarize_continuous(results, version = "complex_covars")
#cont <- summarize_continuous(results, version = "main_terms")
#cont <- summarize_continuous(results, version = "overfit")

#generate_bias_plot(cont)




##########################################
# MSE
##########################################
generate_MSE_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "mse")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("_mse"), names_to = "method", values_to = "MSE") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  plotdat_hi <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("_himse"), names_to = "method", values_to = "himse") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  plotdat_lo <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("_lomse"), names_to = "method", values_to = "lomse") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  pdci <- inner_join(plotdat_lo, plotdat_hi)
  
  if(exclude_clustered){
    pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot(data = plotdat,
         aes(x = ranef_sd,
             group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`)),
         color = method
  ) +
    geom_point(aes(y = MSE,
                   color = method,
                   shape = `Treatment effect`),
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(data = pdci,
                  aes(x = ranef_sd,
                      color = method,
                      group = interaction(`Treatment effect`, method, ranef_sd, `Control arm modeled as`),
                      ymin = lomse, ymax = himse,
                      linetype = `Control arm modeled as`),
                  width = .3, position = position_dodge(width = dodge_width)) +
    facet +
    theme(legend.position = "none",
          title = element_text(size = subtitle_size),
          axis.text = element_text(size = axis_text_size),
          #strip.text = element_text(size = strip_text_size),
          strip.text = ggtext::element_markdown(size = strip_text_size),
          axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept", title = "Mean squared error", y = NULL)
}
#cont <- summarize_continuous(results)
#bin <- summarize_logit(results)
#generate_MSE_plot(cont)
#generate_MSE_plot(bin)





# dwbin <- .6
# dwcont <- .6
# axis_title_size <- 8

##########################################
# TIE
##########################################
generate_TIE_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "TIE")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>% filter(t_eff == 0) %>% 
    pivot_longer(cols = contains("TIE"), names_to = "method", values_to = "TIE") %>%
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           id = row_number()) %>% group_by(id) %>% 
    mutate(TIE_lo = prop.test(x = TIE*nsims, n = nsims)$conf.int[1],
           TIE_hi = prop.test(x = TIE*nsims, n = nsims)$conf.int[2]) %>%
    ungroup() %>% select(-id) %>% 
    mutate(method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  if(exclude_clustered){
    #pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot(data = plotdat,
         aes(x = ranef_sd,
             group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`)),
         color = method
  ) +
    geom_hline(aes(yintercept = .05), linetype = 2, alpha = .5) +
    geom_point(aes(y = TIE,
                   color = method),
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(
      linetype = `Control arm modeled as`,
      color = method,
      ymin = TIE_lo, ymax = TIE_hi),
      width = .2, position = position_dodge(width = dodge_width)) +
    facet +
    ylim(0.01, .10) +
    theme(legend.position = "none",
          title = element_text(size = subtitle_size),
          axis.text = element_text(size = axis_text_size),
          #strip.text = element_text(size = strip_text_size),
          strip.text = ggtext::element_markdown(size = strip_text_size),
          axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept", title = "Type I error", y = NULL)
}

#generate_TIE_plot(cont)


##########################################
# coverage
##########################################
generate_coverage_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "coverage")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>% 
    pivot_longer(cols = contains("_cov"), names_to = "method", values_to = "coverage") %>%
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           id = row_number()) %>% group_by(id) %>% 
    mutate(cov_lo = prop.test(x = coverage*nsims, n = nsims)$conf.int[1],
           cov_hi = prop.test(x = coverage*nsims, n = nsims)$conf.int[2]) %>%
    ungroup() %>% select(-id) %>% 
    mutate(method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  
  if(exclude_clustered){
    #pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot(data = plotdat,
         aes(x = ranef_sd,
             group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`)),
         color = method
  ) +
    geom_hline(aes(yintercept = .95), linetype = 2, alpha = .5) +
    geom_point(aes(y = coverage,
                   color = method,
                   shape = `Treatment effect`),
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(ymin = cov_lo, ymax = cov_hi,
                      color = method,
                      linetype = `Control arm modeled as`),
                  width = .2, position = position_dodge(width = dodge_width)) +
    facet +
    ylim(.9, 1) +
    theme(legend.position = "none",
          title = element_text(size = subtitle_size),
          axis.text = element_text(size = axis_text_size),
          #strip.text = element_text(size = strip_text_size),
          strip.text = ggtext::element_markdown(size = strip_text_size),
          axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept", title = "Coverage", y = NULL)
}
#generate_coverage_plot(cont)

##########################################
# power
##########################################
generate_power_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "power")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>% filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("power"), names_to = "method", values_to = "power") %>%
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           id = row_number()) %>% group_by(id) %>% 
    mutate(power_lo = prop.test(x = power*nsims, n = nsims)$conf.int[1],
           power_hi = prop.test(x = power*nsims, n = nsims)$conf.int[2]) %>%
    ungroup() %>% select(-id) %>% 
    mutate(method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  if(exclude_clustered){
    #pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot(data = plotdat,
         aes(x = ranef_sd,
             group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`)),
         color = method
  ) +
    geom_point(aes(y = power,
                   color = method),
               #shape = `Treatment effect`),
               shape = 17,
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(ymin = power_lo, ymax = power_hi,
                      color = method,
                      linetype = `Control arm modeled as`),
                  width = .3, position = position_dodge(width = dodge_width)) +
    facet +
    theme(legend.position = "none",
          title = element_text(size = subtitle_size),
          axis.text = element_text(size = axis_text_size),
          #strip.text = element_text(size = strip_text_size),
          strip.text = ggtext::element_markdown(size = strip_text_size),
          axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept", title = "Power", y = NULL)
}
#generate_power_plot(cont)
#generate_power_plot(bin)




##########################################
# mean std error of estimate
##########################################
generate_se_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "se")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("mean_se"), names_to = "method", values_to = "mean_SE") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  plotdat_hi <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("_cise_hi"), names_to = "method", values_to = "cise_hi") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  plotdat_lo <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("_cise_lo"), names_to = "method", values_to = "cise_lo") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  pdci <- inner_join(plotdat_lo, plotdat_hi)
  
  if(exclude_clustered){
    pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot() +
    geom_point(data = plotdat,
               aes(x = ranef_sd,
                   y = mean_SE,
                   group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`),
                   color = method,
                   shape = `Treatment effect`),
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(data = pdci,
                  aes(x = ranef_sd,
                      color = method,
                      group = interaction(`Treatment effect`, method, ranef_sd, `Control arm modeled as`),
                      ymin = cise_lo, ymax = cise_hi,
                      linetype = `Control arm modeled as`),
                  width = .3,
                  position = position_dodge(width = dodge_width)) +
    facet +
    theme(legend.position = "none",
          title = element_text(size = subtitle_size),
          axis.text = element_text(size = axis_text_size),
          #strip.text = element_text(size = strip_text_size),
          strip.text = ggtext::element_markdown(size = strip_text_size),
          axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept", title = "Mean estimated SE", y = NULL)
}
#generate_se_plot(cont)



##########################################
# SD of point estimates
##########################################
generate_pe_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "pe_sd")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>%  
    pivot_longer(cols = contains("pe_sd"), names_to = "method", values_to = "sd_point_est") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           id = row_number()) %>% group_by(id) %>% 
    mutate(mean_SD_lo = sd_point_est - 1.96*sd_point_est/sqrt(2*(nsims - 1)),
           mean_SD_hi = sd_point_est + 1.96*sd_point_est/sqrt(2*(nsims - 1))) %>%
    ungroup() %>% select(-id) %>% 
    mutate(method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  if(exclude_clustered){
    #pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot(data = plotdat,
         aes(x = ranef_sd,
             group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`)),
         color = method
  ) +
    geom_point(aes(y = sd_point_est,
                   color = method,
                   shape = `Treatment effect`),
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(x = ranef_sd,
                      color = method,
                      group = interaction(`Treatment effect`, method, ranef_sd, `Control arm modeled as`),
                      ymin = mean_SD_lo, ymax = mean_SD_hi,
                      linetype = `Control arm modeled as`),
                  width = .3, position = position_dodge(width = dodge_width)) +
    
    facet + 
    theme(#legend.position = "none",
      title = element_text(size = subtitle_size),
      axis.text = element_text(size = axis_text_size),
      #strip.text = element_text(size = strip_text_size),
      strip.text = ggtext::element_markdown(size = strip_text_size),
      axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept", title = "Empirical SD of point estimates", y = NULL)
}
#generate_pe_plot(cont)



##########################################
# mean SE divided by SD of point estimates
##########################################
generate_var_ratio_plot <- function(dat, dodge_width = .4, exclude_clustered = F){
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style, type = "ratio")
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  if(dat$link_function[1] %in% c("logit","log")){
    dodge_width <- dwbin
  } else {
    dodge_width <- dwcont
  } 
  
  plotdat <- dat %>% #filter(t_eff != 0) %>% 
    pivot_longer(cols = contains("_var_ratio"), names_to = "method", values_to = "emp_se_div_sd_point_est") %>% 
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  if(exclude_clustered){
    #pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  ggplot(data = plotdat,
         aes(x = ranef_sd,
             group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`),
             #fill = interaction(`Treatment effect`, method),
             #color = interaction(`Treatment effect`, method))
             color = method,
             fill = method)
  ) +
    geom_hline(aes(yintercept = 1), linetype = 2, alpha = .5) +
    geom_point(aes(y = emp_se_div_sd_point_est,
                   #color = method,
                   #fill = interaction(`Treatment effect`, method),
                   #color = interaction(`Treatment effect`, method),
                   shape = interaction(`Treatment effect`,`Control arm modeled as`)),
               position = position_dodge(width = dodge_width)) +
    facet +    
    ylim(0.92, 1.07) +
    scale_shape_manual(values = c(16,17, 16, 17)) +
    theme(legend.position = "none",
          title = element_text(size = subtitle_size),
          axis.text = element_text(size = axis_text_size),
          #strip.text = element_text(size = strip_text_size),
          strip.text = ggtext::element_markdown(size = strip_text_size),
          axis.title = element_text(size = axis_title_size)) +
    labs(x = "Standard deviation of treatment arm random intercept",
         title = "Ratio of estimated SE to empirical SD", y = NULL)
}
#generate_var_ratio_plot(cont)


##########################################
# empty plot for last panel
##########################################
generate_legend_plot <- function(dat, exclude_clustered = F){
  
  style <- dat[["style"]][[1]]
  
  facet_plan <- get_facet(dat = dat, style = style)
  dat <- facet_plan$dat
  facet <- facet_plan$facet
  
  plotdat <- dat %>%
    pivot_longer(cols = contains("mean_se"), names_to = "method", values_to = "power") %>%
    mutate(`Control arm modeled as` = ifelse(str_detect(method, "_cl"),"clustered", "unclustered"),
           method = substr(method, start = 1, stop = 4),
           method = ifelse(method == "gee_", "GEE", ifelse(method == "glmm", "GLMM", "TMLE")))
  
  if(exclude_clustered){
    #pdci <- pdci %>% filter(`Control arm modeled as` == "unclustered")
    plotdat <- plotdat %>% filter(`Control arm modeled as` == "unclustered")
  }
  
  
  cowplot::get_legend(ggplot(data = plotdat,
                             aes(x = ranef_sd#,
                                 #group = interaction(`Treatment effect`, ranef_sd, method, `Control arm modeled as`)
                             ),
                             color = method) +
                        geom_point(aes(y = power,
                                       color = method,
                                       shape = `Treatment effect`)) +
                        geom_errorbar(aes(ymin = 0, ymax = 1#,
                                          #linetype = `Control arm modeled as`
                        )) +
                        facet +
                        theme(legend.box = "horizontal",
                              legend.key.size = unit(legend_icon_size,"line"),
                              legend.text = element_text(size = legend_size),
                              legend.title = element_text(size = legend_size)) + 
                        guides(colour = guide_legend(override.aes = list(size = legend_dot_size))) +
                        guides(shape = guide_legend(override.aes = list(size = legend_dot_size)))
  )
}




















