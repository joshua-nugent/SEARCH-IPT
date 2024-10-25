source("2_run_sims.R")


sizes_irgtt <- list(
  c(cluster_size_t = 5, clusters_per_arm_t = 12,
    cluster_size_c = 5, clusters_per_arm_c = 12,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 5, clusters_per_arm_t = 30,
    cluster_size_c = 5, clusters_per_arm_c = 30,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 5, clusters_per_arm_t = 50,
    cluster_size_c = 5, clusters_per_arm_c = 50,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 10, clusters_per_arm_t = 12,
    cluster_size_c = 10, clusters_per_arm_c = 12,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 10, clusters_per_arm_t = 30,
    cluster_size_c = 10, clusters_per_arm_c = 30,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 10, clusters_per_arm_t = 50,
    cluster_size_c = 10, clusters_per_arm_c = 50,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 20, clusters_per_arm_t = 12,
    cluster_size_c = 20, clusters_per_arm_c = 12,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 20, clusters_per_arm_t = 30,
    cluster_size_c = 20, clusters_per_arm_c = 30,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 20, clusters_per_arm_t = 50,
    cluster_size_c = 20, clusters_per_arm_c = 50,
    control_arm_irgtt = T, control_arm_no_ranef = T))



nsims <- 2000
core_n <- 9

Sys.time()
results_irgtt <- run_full_set(cores = core_n, seed = 11,
                              sizes = sizes_irgtt,
                              link_function = "identity",
                              nsims = nsims)
saveRDS(results_irgtt, "results/results_irgtt.rds")
Sys.time()
rm(results_irgtt)
gc()

Sys.time()
results_irgtt_binary <- run_full_set(cores = core_n, seed = 11,
                                     sizes = sizes_irgtt,
                                     t_eff = c(.35, 0), # make larger for binary
                                     link_function = "logit",
                                     nsims = nsims)
saveRDS(results_irgtt_binary, "results/results_irgtt_binary.rds")
Sys.time()
rm(results_irgtt_binary)
gc()


Sys.time()
results_irgtt_misspec <- run_full_set(cores = core_n, seed = 11,
                                      tmle_misspec_test = T,
                                      sizes = sizes_irgtt,
                                      link_function = "identity",
                                      nsims = nsims)
saveRDS(results_irgtt_misspec, "results/results_irgtt_misspec.rds")
Sys.time()
rm(results_irgtt_misspec)
gc()





sizes_unbalanced <- list(
  # .5x smaller N in control arm
  c(cluster_size_t = 5, clusters_per_arm_t = 30,
    cluster_size_c = 5, clusters_per_arm_c = 15,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 5, clusters_per_arm_t = 50,
    cluster_size_c = 5, clusters_per_arm_c = 25,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 10, clusters_per_arm_t = 30,
    cluster_size_c = 10, clusters_per_arm_c = 15,
    control_arm_irgtt = T, control_arm_no_ranef = T),
  c(cluster_size_t = 10, clusters_per_arm_t = 50,
    cluster_size_c = 10, clusters_per_arm_c = 25,
    control_arm_irgtt = T, control_arm_no_ranef = T))


nsims <- 2000
core_n <- 9


Sys.time()
results_irgtt_unbalanced <- run_full_set(cores = core_n, seed = 11,
                              sizes = sizes_unbalanced,
                              link_function = "identity",
                              nsims = nsims)
saveRDS(results_irgtt_unbalanced, "results/results_irgtt_unbalanced.rds")
Sys.time()
rm(results_irgtt_unbalanced)
gc()

Sys.time()
results_irgtt_binary_unbalanced <- run_full_set(cores = core_n, seed = 11,
                                     sizes = sizes_unbalanced,
                                     t_eff = c(.35, 0), # make larger for binary
                                     link_function = "logit",
                                     nsims = nsims)
saveRDS(results_irgtt_binary_unbalanced, "results/results_irgtt_binary_unbalanced.rds")
Sys.time()
rm(results_irgtt_binary_unbalanced)
gc()
plan(sequential)

















