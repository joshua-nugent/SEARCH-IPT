library(tidyverse)
library(ltmle)
library(future)
library(doFuture)
library(foreach)
library(sandwich)
library(lme4)
library(gee)
source("0_DGP_and_helpers.R")
source("1_fit_models.R")

one_set <- function(specification, nsims = 4, cores = 8, trial = T, seed = NA,
                    tmle_misspec_test = F,
                    estimate_pscore = T){
  
  sizes <- specification[[1,"sizes"]]
  sizes <- sizes[[1]]
  #print(sizes)
  if(!is.na(seed)){
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }
  
  if(specification[[1,"link_function"]] == "identity"){
    PATE <- cbind.data.frame(target_parameter = specification[[1,"t_eff"]])
  } else { # Get true PATE values by averaging...
    pop <- bind_rows(lapply(as.list(1:500), function(x) gen_obs_dat(LPS_only = T,
                                                                     trial = trial,
                                                                     sizes = sizes,#specification[[1,"sizes"]],
                                                                     t_eff = specification[[1,"t_eff"]],
                                                                     #control_arm_no_ranef = specification[[1,"control_arm_no_ranef"]],
                                                                     ranef_sd = specification[[1,"ranef_sd"]],
                                                                     link_function = specification[[1,"link_function"]],
                                                                     sim_type = specification[[1,"sim_type"]])))
    UY <- runif(nrow(pop))
    P0 <- mean(as.numeric(plogis(pop$LP0) > UY))
    P1 <- mean(as.numeric(plogis(pop$LP1) > UY))
    PATE <- cbind.data.frame(target_parameter = (P1/(1-P1)) / (P0/(1-P0)))
  }
  
  #print(sizes)
  #print(specification)
  plan(multisession, workers = cores)
  if(!is.na(seed)){
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }
  results <- data.frame(foreach(j = 1:nsims, .combine = rbind,
                                .options.future = list(seed = T, globals = structure(T)),
                                .errorhandling = "remove") %dofuture% {
    fit_models(dat = gen_obs_dat(sizes = sizes,
                                 trial = trial,
                                 #control_arm_no_ranef = sizes[[1,"control_arm_no_ranef"]],
                                 t_eff = specification[[1,"t_eff"]],
                                 ranef_sd = specification[[1,"ranef_sd"]],
                                 link_function = specification[[1,"link_function"]],
                                 sim_type = specification[[1,"sim_type"]]),
               sizes = sizes,
               tmle_misspec_test = tmle_misspec_test,
               #control_arm_irgtt = sizes[[1]][["control_arm_irgtt"]],
               link_function = specification[[1,"link_function"]],
               #v = specification[[1,"v"]],
               sim_type = specification[[1,"sim_type"]],
               estimate_pscore = estimate_pscore)
  })
  plan(sequential)
  return(cbind.data.frame(nsims = nsims,
                          specification,
                          control_arm_irgtt = sizes[["control_arm_irgtt"]],
                          control_arm_no_ranef = sizes[["control_arm_no_ranef"]],
                          control_arm_clusters = sizes[["clusters_per_arm_c"]],
                          treatment_arm_clusters = sizes[["clusters_per_arm_t"]],
                          control_arm_cluster_size = sizes[["cluster_size_c"]],
                          treatment_arm_cluster_size = sizes[["cluster_size_t"]],
                          results, PATE))
}


run_full_set <- function(sim_type = c("complex_covars", "main_terms", "overfit"),
                         link_function = "identity",
                         t_eff = c(.25, 0),
                         ranef_sd = c(icc_to_sd(.001),
                                      icc_to_sd(.01),
                                      icc_to_sd(.05),
                                      icc_to_sd(.10)),
                         #v = 10,
                         tmle_misspec_test = F,
                         estimate_pscore = T, nsims = 4, cores = 8, seed = NA, trial = T,
                         sizes = list(c(cluster_size_t = 2, clusters_per_arm_t = 10,
                                        cluster_size_c = 2, clusters_per_arm_c = 10,
                                        control_arm_irgtt = T, control_arm_no_ranef = T))){

  specification <- expand_grid(sizes,
                               sim_type,
                               link_function, 
                               t_eff, 
                               ranef_sd,
                               #v,
                               seed)
  #print(specification)
  #return(specification)
  
  return(data.frame(foreach(j = 1:nrow(specification), .combine = rbind) %do%
                      {one_set(specification = specification[j,],
                               seed = seed,
                               trial = trial,
                               cores = cores,
                               nsims = nsims,
                               tmle_misspec_test = tmle_misspec_test,
                               estimate_pscore = estimate_pscore)}))
}
# run_full_set(link_function = "logit", tmle_misspec_test = T)
# sim_type = "nonadditive"
# link_function = "identity"
# t_eff = c(.25, 0)
# ranef_sd = c(.05)
# control_arm_no_ranef = T
# v = 10
# control_arm_irgtt = T
# estimate_pscore = T
# nsims = 2
# cores = 10
# seed = NA
# 
# sizes <- list(c(cluster_size_t = 5, clusters_per_arm_t = 10,
#                 cluster_size_c = 5, clusters_per_arm_c = 10,
#                 control_arm_irgtt = T, control_arm_no_ranef = T))
# 
# specification <- expand_grid(sizes,
#                              sim_type,
#                              link_function,
#                              t_eff,
#                              ranef_sd,
#                              control_arm_no_ranef,
#                              v,
#                              control_arm_irgtt,
#                              seed)
# specification[1,]
# 
# one_set(specification = specification[1,],
#         cores = cores,
#         nsims = nsims,
#         estimate_pscore = estimate_pscore)
# 
# 
# fit_models(dat = gen_obs_dat(seed = specification[[1,"seed"]],
#                              sizes = specification[[1,"sizes"]],
#                              control_arm_no_ranef = specification[[1,"control_arm_no_ranef"]],
#                              t_eff = specification[[1,"t_eff"]],
#                              ranef_sd = specification[[1,"ranef_sd"]],
#                              link_function = specification[[1,"link_function"]],
#                              sim_type = specification[[1,"sim_type"]]),
#            control_arm_irgtt = specification[[1,"control_arm_irgtt"]],
#            link_function = specification[[1,"link_function"]],
#            v = specification[[1,"v"]],
#            estimate_pscore = T)
