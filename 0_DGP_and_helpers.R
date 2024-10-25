require(tidyverse)
`%nin%` <- negate(`%in%`)

icc_to_sd <- function(icc = .01){
  sqrt(icc / (1-icc))
}

bin_icc <- function(s2 = .05^2){
  s2 / (s2 + (pi^2/3))
}

# see https://academic.oup.com/ije/article/52/4/1276/7076266
get_folds <- function(n_eff){
  if(n_eff < 30){
    return(n_eff)
  } else if(n_eff < 500){
    return(20)
  } else if(n_eff < 5000){
    return(10)
  } else if(n_eff < 10000){
    return(5)
  } else {
    return(4)
  }
}


# see https://academic.oup.com/ije/article/52/4/1276/7076266
get_n_eff <- function(outcome_type, variable_vector){
  n <- length(variable_vector)
  if(outcome_type == "binomial"){
    py1 <- sum(variable_vector)
    py0 <- n - py1
    nrare <- min(py1, py0)
    return(min(n, 5*nrare))
  } else {
    return(n)
  }
}

gen_obs_dat <- function(link_function = "identity",
                        trial = T,
                        sim_type = "complex_covars",
                        sizes = list(c(cluster_size_t = 3, clusters_per_arm_t = 10,
                                       cluster_size_c = 3, clusters_per_arm_c = 10,
                                       control_arm_no_ranef = T)),
                        t_eff = .25, ranef_sd = .1, seed = NA,
                        LPS_only = F){
  
  if(length(sizes) == 1){
    sizes <- sizes[[1]]
  }
  
  control_arm_no_ranef <- sizes[["control_arm_no_ranef"]] # is U_E for control arm 0 or not?
  
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  # Set up units
  total_units_t <- sizes[["cluster_size_t"]] * sizes[["clusters_per_arm_t"]]
  total_units_c <- sizes[["cluster_size_c"]] * sizes[["clusters_per_arm_c"]]
  total_units <- total_units_t + total_units_c
  clust_id_t <- rep(1:sizes[["clusters_per_arm_t"]], each = sizes[["cluster_size_t"]])
  clust_id_c <- rep((max(clust_id_t) + 1):(max(clust_id_t) + sizes[["clusters_per_arm_c"]]), each = sizes[["cluster_size_c"]])
  clust_id_c_irgtt <- (max(clust_id_t) + 1):(max(clust_id_t) + total_units_c)
  clust_id <- c(clust_id_t, clust_id_c)
  clust_id_irgtt <- c(clust_id_t, clust_id_c_irgtt)

  X_bin <- rbinom(n = total_units, size = 1, prob = .25)
  if(trial){
    A <- c(rep(1, times = total_units_t), rep(0, times = total_units_c))
    X_cont <- rnorm(n = total_units, mean = 0, sd = 1)
  } else {
    A <- c(rep(1, times = total_units_t), rep(0, times = total_units_c))
    X_cont <- rnorm(n = total_units, mean = .1*A, sd = 1)
  }

  if(ranef_sd == 0){
    UE <- 0
  } else {
    if(control_arm_no_ranef){
      UEs_t <- rnorm(n = sizes[["clusters_per_arm_t"]], mean = 0, sd = ranef_sd)
      UE <- c(rep(UEs_t, each = sizes[["cluster_size_t"]]),
              rep(0, times = total_units_c))
    } else {
      UEs_t <- rnorm(n = sizes[["clusters_per_arm_t"]], mean = 0, sd = ranef_sd)
      UE_t <- c(rep(UEs_t, each = sizes[["cluster_size_t"]]))
      UEs_c <- rnorm(n = sizes[["clusters_per_arm_c"]], mean = 0, sd = ranef_sd)
      UE_c <- c(rep(UEs_c, each = sizes[["cluster_size_c"]]))
      UE <- c(UE_t, UE_c)
      clust_id_irgtt <- "none"
    }
  }
  
  
  if(link_function == "logit"){
    intercept <- -5
    
    if(sim_type == "main_terms"){
      LP0 <- 4.5 + intercept + .5*X_bin + X_cont + UE
      LP1 <- 4.5 + intercept + .5*X_bin + X_cont + UE + t_eff
    } else if(sim_type == "complex_covars"){
      LP0 <- 3 + intercept - .5*X_bin + X_cont^2 + (X_cont^2)*X_bin + UE
      LP1 <- 3 + intercept - .5*X_bin + X_cont^2 + (X_cont^2)*X_bin + UE + t_eff
    } else if(sim_type == "overfit"){
      LP0 <- 5 + intercept + UE
      LP1 <- 5 + intercept + UE + t_eff
    } else if(sim_type == "nonadditive"){
      
      
    } else {
      print("Simulation structure not entered correctly. Must be main_terms, complex_covars, or overfit")
    }
    LP0_marg <- LP0 - UE
    LP1_marg <- LP1 - UE
    P0 <- plogis(LP0)
    P1 <- plogis(LP1)

    UY <- runif(total_units)
    Y0 <- as.numeric(P0 > UY)
    Y1 <- as.numeric(P1 > UY)
    Y <- ifelse(A == 1, Y1, Y0)
    
  } else {
    # Continuous outcomes
    UY <- rnorm(n = total_units, sd = 1)
    
    if(sim_type == "main_terms"){
      LP0 <- .5*X_bin + X_cont + UY + UE
      LP1 <- .5*X_bin + X_cont + UY + UE + t_eff
    } else if(sim_type == "complex_covars"){
      LP0 <- -.5*X_bin + X_cont^2 + (X_cont^2)*X_bin + UY + UE
      LP1 <- -.5*X_bin + X_cont^2 + (X_cont^2)*X_bin + UY + UE + t_eff
    } else if(sim_type == "overfit"){
      LP0 <- UY + UE
      LP1 <- UY + UE + t_eff
    } else if(sim_type == "nonadditive"){

    } else {
      print("Simulation structure not entered correctly. Must be main_terms, complex_covars, overfit, or non-additive")
    }
    
    Y <- ifelse(A == 1, LP1, LP0)  
  }
  
  if(LPS_only){
    return(cbind.data.frame(LP0 = LP0, LP1 = LP1, t_eff))
  }
  return(cbind.data.frame(UE, clust_id, clust_id_irgtt, A, X_cont, X_bin, Y))
}
# head(gen_obs_dat(sizes = list(c(cluster_size_t = 2, clusters_per_arm_t = 10,
#                            cluster_size_c = 2, clusters_per_arm_c = 10,
#                            control_arm_irgtt = F,
#                            control_arm_no_ranef = F)),
#             link_function = "identity", sim_type = "nonadditive"))




