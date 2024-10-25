library(tidyverse)
library(ltmle)
library(SuperLearner)
library(sandwich)
library(lme4)
library(gee)
library(lmerTest)
library(lmtest)

HouseholdIC <- function(recordIC, id) {
  if (is.null(id)) return(recordIC)
  householdIC <- as.matrix(aggregate(recordIC, list(id=id), sum)[, -1, drop=FALSE])
  num.records <- nrow(recordIC)
  num.households <- nrow(householdIC)
  householdIC <- householdIC * num.households / num.records
  return(householdIC)
}

fit_tmle <- function(dat,
                     tmle_id_irgtt = T,
                     link_function = "identity",
                     verbose = T,
                     sizes,
                     tmle_break_CV = F,
                     tmle_ignore_clustering = F,
                     estimate_pscore = T,
                     sim_type = "nonadditive",
                     SLL = c("SL.mean", "SL.glm","SL.earth")){
  n <- nrow(dat)
  lower_g_bound <- 5/sqrt(n)/log(n) # following TMLE recommendation
  
  if(link_function == "identity"){
    n_eff <- get_n_eff(outcome_type = "continuous", variable_vector = dat[["Y"]])
  } else if(link_function == "logit"){
    n_eff <- get_n_eff(outcome_type = "binomial", variable_vector = dat[["Y"]])
  }
  SL_cv_folds <- get_folds(n_eff)
  
  if(estimate_pscore){
    gform <- NULL
  } else {
    gform <- "A ~ 1"
  }
  
  if(tmle_id_irgtt){
    tmle_ids <- dat$clust_id_irgtt
  } else {
    tmle_ids <- dat$clust_id
  }
  
  nclust <- length(unique(dat$clust_id_irgtt))
  nclust_cl <- length(unique(dat$clust_id))
  
  if(tmle_ignore_clustering){
    if(verbose){print("generating estimates for TMLE ignoring clustering. Bad!")}
    tmle_ids <- NULL
    nclust <- n
  }
  
  
  if(tmle_break_CV){
    if(verbose){print("generating estimates for continuous outcome (TMLE) without proper CV structure. Bad!")}
    
    ltmle_mod_0 <- suppressWarnings(suppressMessages(
      ltmle(data = dat %>% dplyr::select(all_of(c("X_cont", "X_bin", "A", "Y"))),
            Ynodes = "Y", Anodes = "A",
            gbounds = c(lower_g_bound, 1),
            variance.method = "ic",
            SL.library = SLL,
            SL.cvControl = list(V = SL_cv_folds),
            abar = 0)))
    ltmle_mod_1 <- suppressWarnings(suppressMessages(
      ltmle(data = dat %>% dplyr::select(all_of(c("X_cont", "X_bin", "A", "Y"))),
            Ynodes = "Y", Anodes = "A",
            gbounds = c(lower_g_bound, 1),
            variance.method = "ic",
            SL.library = SLL,
            SL.cvControl = list(V = SL_cv_folds),
            abar = 1)))
    
    ltmle_pe <- ltmle_mod_1$estimates[["tmle"]] - ltmle_mod_0$estimates[["tmle"]]
    ate_ic <- ltmle_mod_1$IC$tmle - ltmle_mod_0$IC$tmle
    
    
    householdIC <- aggregate(ate_ic, list(id=tmle_ids), sum)[, -1, drop=FALSE]
    num.records <- length(ate_ic)
    num.households <- nrow(householdIC)
    householdIC <- unlist(householdIC * num.households / num.records)
    
    ltmle_se <- sqrt(var(householdIC)/nclust)
    ltmle_lo <- ltmle_pe - ltmle_se*qt(.975, df = nclust - 2)
    ltmle_hi <- ltmle_pe + ltmle_se*qt(.975, df = nclust - 2)
    
  } else {
    if(link_function == "logit") {
      if(verbose){print("generating estimates for binary outcome (TMLE)")}
      
      ltmle_mod <- suppressWarnings(suppressMessages(
        ltmle(data = dat %>% dplyr::select(all_of(c("X_cont", "X_bin", "A", "Y"))),
              Ynodes = "Y", Anodes = "A",
              id = tmle_ids,
              gbounds = c(lower_g_bound, 1),
              variance.method = "ic", 
              SL.library = SLL,
              SL.cvControl = list(V = SL_cv_folds),
              abar = list(1,0))))
      summary_ltmle_mod <- summary(ltmle_mod)
      ltmle_pe <- summary_ltmle_mod$effect.measures$OR$estimate
      ltmle_se <- summary_ltmle_mod$effect.measures$OR$std.dev
      ltmle_lo <- exp(log(ltmle_pe) - ltmle_se*qt(.975, df = nclust - 2))
      ltmle_hi <- exp(log(ltmle_pe) + ltmle_se*qt(.975, df = nclust - 2))
    } else {
      if(verbose){print("generating estimates for continuous outcome (TMLE)")}
      
      ltmle_mod <- suppressWarnings(suppressMessages(
        ltmle(data = dat %>% dplyr::select(all_of(c("X_cont", "X_bin", "A", "Y"))),
              Ynodes = "Y", Anodes = "A",
              id = tmle_ids,
              gbounds = c(lower_g_bound, 1),
              variance.method = "ic",
              SL.library = SLL,
              SL.cvControl = list(V = SL_cv_folds),
              abar = list(1,0))))
      summary_ltmle_mod <- summary(ltmle_mod)
      ltmle_pe <- summary_ltmle_mod$effect.measures$ATE$estimate
      ltmle_se <- summary_ltmle_mod$effect.measures$ATE$std.dev
      ltmle_lo <- ltmle_pe - ltmle_se*qt(.975, df = nclust - 2)
      ltmle_hi <- ltmle_pe + ltmle_se*qt(.975, df = nclust - 2)
    }
  }
  output <- cbind.data.frame(ltmle_pe = ltmle_pe, ltmle_se = ltmle_se,
                             ltmle_lo = ltmle_lo, ltmle_hi = ltmle_hi)
  if(tmle_ignore_clustering){
    colnames(output) <- paste0(colnames(output), "_ignore_clustering")
  }
  if(tmle_break_CV){
    colnames(output) <- paste0(colnames(output), "_break_CV")
  }
  return(output)
}

fit_gee <- function(dat,
                    link_function = "identity",
                    gee_ind = F,
                    verbose = T,
                    sizes){
  library(gee)
  library(sandwich)
  
  if(link_function == "logit") {
    if(verbose){print("generating estimates for logit GEE")}
    if(gee_ind){
      gee <- tryCatch(
        glm(Y ~ A + X_cont + X_bin,
            data = dat, family = binomial(link = "logit")) %>% coeftest(vcov = sandwich),
        error = function(e){
          print(e$message)
          return(NULL)
        }  
      )
    } else {
      gee <- tryCatch(
        gee(Y ~ A + X_cont + X_bin, id = dat$clust_id_irgtt, corstr = "exchangeable",
            maxiter = 50,
            data = dat, family = binomial(link = "logit")),
        error = function(e){
          print(e$message)
          return(NULL)
        }
      )
      
      if(is.null(gee)){
        gee_pe <- NA
        gee_se <- NA
        gee_lo <- NA
        gee_hi <- NA
      } else {
        if(gee_ind){
          gee_pe <- exp(gee["A","Estimate"])
          gee_se <- gee["A","Std. Error"]
          gee_lo <- exp(log(gee_pe) - qnorm(0.975) * gee_se)
          gee_hi <- exp(log(gee_pe) + qnorm(0.975) * gee_se)
        } else {
          gee_pe <- exp(summary(gee)$coefficients["A","Estimate"])
          gee_se <- summary(gee)$coefficients["A","Robust S.E."]
          gee_lo <- exp(log(gee_pe) - qnorm(0.975) * gee_se)
          gee_hi <- exp(log(gee_pe) + qnorm(0.975) * gee_se)
        }
      }
    }
  } else {
    if(verbose){print("generating estimates for continuous GEE")}
    gee <- gee(Y ~ A + X_cont + X_bin, id = dat$clust_id_irgtt, corstr = "exchangeable",
               data = dat)
    gee_pe <- summary(gee)$coefficients["A","Estimate"]
    gee_se <- summary(gee)$coefficients["A","Robust S.E."]
    gee_lo <- gee_pe - qnorm(0.975) * gee_se
    gee_hi <- gee_pe + qnorm(0.975) * gee_se
  }
  if(verbose){print("GEE fitting complete")}
  return(cbind.data.frame(gee_pe = gee_pe,
                          gee_se = gee_se,
                          gee_lo = gee_lo,
                          gee_hi = gee_hi))
}

fit_glmm <- function(dat,
                     link_function = "identity",
                     verbose = T,
                     sizes){
  library(lme4)
  library(lmerTest)
  
  if(link_function == "logit") {
    glmm <- tryCatch(
      glmer(Y ~ A + X_cont + X_bin +
              (0 + A | clust_id_irgtt),
            control = glmerControl(optCtrl = list(maxit = 1e8, maxfun = 1e8),
                                   boundary.tol = 1e-6),
            family = binomial(link = "logit"),
            data = dat),
      error = function(e){
        print(e$message)
        return(NULL)
      }
    )
    if(is.null(glmm)){
      glmm_pe <- NA
      glmm_se <- NA
      glmm_lo <- NA
      glmm_hi <- NA
    } else {
      glmm_pe <- exp(summary(glmm)$coefficients["A","Estimate"])
      glmm_se <- summary(glmm)$coefficients["A","Std. Error"]
      glmm_lo <- exp(log(glmm_pe) - qnorm(0.975) * glmm_se)
      glmm_hi <- exp(log(glmm_pe) + qnorm(0.975) * glmm_se)
    }
  } else {
    # Note using Nugent paper to find optimal DF for LMM in this situation (Satterthwaite)
    glmm <- lmer(Y ~ A + X_cont + X_bin + (0 + A | clust_id_irgtt),
                 data = dat, REML = T)
    glmm_pe <- summary(glmm)$coefficients["A","Estimate"]
    glmm_se <- summary(glmm)$coefficients["A","Std. Error"]
    satterthwaite_df <- summary(glmm)$coefficients["A","df"]
    glmm_lo <- glmm_pe - qt(0.975, df = satterthwaite_df) * glmm_se
    glmm_hi <- glmm_pe + qt(0.975, df = satterthwaite_df) * glmm_se
  }
  if(verbose){print("GLMM fitting complete")}
  return(cbind.data.frame(glmm_pe = glmm_pe,
                          glmm_se = glmm_se,
                          glmm_lo = glmm_lo,
                          glmm_hi = glmm_hi))
}

fit_models <- function(dat,
                       gee_ind = F,
                       tmle_id_irgtt = T,
                       link_function = "identity",
                       verbose = F,
                       sizes,
                       tmle_misspec_test = F,
                       estimate_pscore = T,
                       sim_type = "nonadditive",
                       SLL = c("SL.mean", "SL.glm","SL.earth")){
  
  tmle_results <- fit_tmle(dat = dat,
                           tmle_id_irgtt = tmle_id_irgtt,
                           link_function = link_function,
                           verbose = verbose,
                           sizes = sizes,
                           estimate_pscore = estimate_pscore,
                           sim_type = sim_type,
                           SLL = SLL)
  
  if(tmle_misspec_test){
    tmle_bad_CV_results <- fit_tmle(dat = dat,
                                    tmle_id_irgtt = tmle_id_irgtt,
                                    link_function = link_function,
                                    verbose = verbose,
                                    sizes = sizes,
                                    tmle_break_CV = T,
                                    estimate_pscore = estimate_pscore,
                                    sim_type = sim_type,
                                    SLL = SLL)
    tmle_bad_IC_results <- fit_tmle(dat = dat,
                                    tmle_id_irgtt = tmle_id_irgtt,
                                    link_function = link_function,
                                    verbose = verbose,
                                    sizes = sizes,
                                    tmle_ignore_clustering = T,
                                    estimate_pscore = estimate_pscore,
                                    sim_type = sim_type,
                                    SLL = SLL)
    return(cbind.data.frame(tmle_results, tmle_bad_CV_results, tmle_bad_IC_results))
  } else {
    gee_results <- fit_gee(dat = dat,
                           link_function = link_function,
                           gee_ind = gee_ind,
                           verbose = verbose,
                           sizes = sizes)
    
    glmm_results <- fit_glmm(dat = dat,
                             link_function = link_function,
                             verbose = verbose,
                             sizes = sizes)
    return(cbind.data.frame(tmle_results, gee_results, glmm_results))
  } 
}


# dat <- gen_obs_dat(sizes = list(c(cluster_size_t = 4, clusters_per_arm_t = 20,
#                                   cluster_size_c = 4, clusters_per_arm_c = 20,
#                                   control_arm_irgtt = T,
#                                   control_arm_no_ranef = T)),
#                    link_function = "identity")
# 
# fit_models(dat = dat,
#            link_function = "identity",
#            verbose = T,
#            tmle_misspec_test = T,
#            sizes = c(cluster_size_t = 4, clusters_per_arm_t = 20,
#                      cluster_size_c = 4, clusters_per_arm_c = 20,
#                      control_arm_irgtt = T,
#                      control_arm_no_ranef = F))




