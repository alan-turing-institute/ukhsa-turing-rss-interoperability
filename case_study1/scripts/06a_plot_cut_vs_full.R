
### Estimate LTLA prevalence ###
library(dplyr)
library(prevdebiasr)
library(parallel)
library(foreach)
library(truncnorm)
source("scripts/SIR_utils.R")

# trans_mats <- readRDS("transmats/poisson_SIR_epi_gamma_1.RDS")
vax_df <- readr::read_csv(file.path(data_dir_to_use, "vaccination.csv"))
pcr_infectious_df <- readr::read_csv(file.path(data_dir_to_use, "moment_match_infectious.csv"))
region_df <- readr::read_csv(file.path(data_dir_to_use, "region.csv")) %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week"))
ltla_df <- readr::read_csv(file.path(data_dir_to_use, "ltla.csv")) %>%
  left_join(vax_df, by = c("ltla", "mid_week")) %>%
  left_join(select(region_df, -c(Nt, nt, M)),
            by = c("phe_region", "mid_week")
  ) %>%
  select(ltla, phe_region, mid_week, Nt, nt, M, V, alpha, beta)


region_df[region_df$phe_region == "London" & region_df$mid_week == "2021-01-17", ]
ltla_df[ltla_df$phe_region == "London" & ltla_df$mid_week == "2021-01-17", ]

alpha_testing <- 3e-4
control_debias <- prevdebiasr::get_control_parameters(
  alpha_testing = alpha_testing
)

reg_df_curr <- region_df[region_df$phe_region == "South West", ]
reg_df_curr[nrow(reg_df_curr) - 10:0, ]
mid_week_unique <- sort(unique(ltla_df$mid_week))

n_cores <- 12
run_type <- c("fast", "full")[1]
clust <- makeCluster(n_cores)
doParallel::registerDoParallel(clust)

out_dir <- paste0(
  "output/AR", control_debias$delta_AR_rho,
  "sd", control_debias$delta_AR_sd,
  "Rsd", control_SIR$R_AR_sd
)
control <- c(control_debias, control_SIR)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(control, file.path(out_dir, "control.RDS"), version = 2)

clusterExport(clust, c("control_debias", "control_SIR"))

### PCR positivity, perfect testing ###
imperfect <- FALSE
type <- "PCR_positive"
region_curr <- c("London", "South West", "South East", "North West")[1]
mid_week_use <- "2021-01-17"
test_df_full <- region_df[region_df$phe_region == region_curr, ]
data_in <- c("single_week", "whole_round")[1]
if(data_in == "single_week") {
  test_df <- test_df_full[test_df_full$mid_week == mid_week_use, ]
}
if(data_in == "whole_round") {
  react_round_8_weeks <- mid_week_unique[mid_week_unique > "2021-01-06" & mid_week_unique < "2021-01-21"]
  test_df_curr_round <- test_df[test_df$mid_week %in% react_round_8_weeks, ]
  test_df_round_aggregated <- test_df_curr_round[1, c("phe_region", "mid_week", "Nt", "nt", "Nr", "nr", "M")]
  test_df_round_aggregated[, c("Nt", "nt", "Nr", "nr")] <-  test_df_curr_round %>% 
    summarise(across(all_of(c("Nt", "nt", "Nr", "nr")), sum))
  test_df_round_aggregated$round <- "8"
  test_df <- test_df_round_aggregated
}


par(mfrow = c(2, 1))
regc <- "Barking and Dagenham"
ltla_df_curr <- as.data.frame(ltla_df[ltla_df$ltla == regc, ])
matplot(1:nrow(ltla_df_curr), cbind(ltla_df_curr$Nt, ltla_df_curr$nt), xaxt = "n", xlab = "")
axis(side = 1, at = 1:nrow(ltla_df_curr), labels = ltla_df_curr$mid_week, las = 2, cex.axis = .7, xlab = "")
matplot(1:nrow(ltla_df_curr), ltla_df_curr$nt / ltla_df_curr$Nt, xaxt = "n", xlab = "")
axis(side = 1, at = 1:nrow(ltla_df_curr), labels = ltla_df_curr$mid_week, las = 2, cex.axis = .7, xlab = "")

j <- 1
M_curr <- test_df$M[j]
control <- control_debias
cut = T
test_df$delta_post_sd <- test_df$delta_post_mean <- NA
nu_approx <- boot::logit((test_df$Nt - test_df$nt)/test_df$M)
n_del <- 2000
n_pi_coarse <- 2000
n_pi_fine_plot <- 2000
pi_range <- c(0, 0.05)
log_bound <- -5000
nu_fixed_pt_convergence_tol <- 1e-3
I_fine <- floor(seq(floor(pi_range[1] * M_curr), floor(pi_range[2] * M_curr), length.out = n_pi_fine_plot))
control$priors$trunc_gauss$upper_trunc_proportion <- .06
extra_I_band <- round(M_curr * .005)
del_seq <- seq(1, 6, len = n_del)
pi_seq_coarse <- seq(0, .05, len = n_pi_coarse)
I_seq_coarse <- unique(floor(M_curr * pi_seq_coarse))

#####################################################
# Evaluate some distributions and define some functions
I_logprior <- log(prevdebiasr:::prior_prevalence(test_df, control))
I_loglik_just_rand <- stats::dhyper(test_df$nr[j], control$I_seq, 
                                    pmax(0, M_curr - control$I_seq), test_df$Nr[j], log = TRUE)
I_logpost <- I_logprior + I_loglik_just_rand
I_post_norm <- exp(I_logpost - max(I_logpost))
I_post_norm <- I_post_norm / sum(I_post_norm)
prior_on_nu_plus_delta <- function(x) exp(-x) * (1 + exp(-x))^(-2)

log_joint_posterior <- function(delta, I, del_seq, nu) {
  approxfun(del_seq + nu[j], prior_on_nu_plus_delta(del_seq + nu[j]))(delta + nu[j]) +
    approxfun(control$I_seq, I_logprior[j, ])(I) +
    approxfun(control$I_seq, I_loglik_just_rand)(I) +
    dbinom(test_df$nt[j], size = I, prob = boot::inv.logit(delta + nu[j]), log = T) +
    dbinom(test_df$Nt[j] - test_df$nt[j], size = M_curr - I, prob = boot::inv.logit(nu[j]), log = T)
}

log_delta_cond_posterior <- function(delta, I, del_seq, nu) {
  approxfun(del_seq + nu[j], prior_on_nu_plus_delta(del_seq + nu[j]))(delta + nu[j]) +
    dbinom(test_df$nt[j], size = I, prob = boot::inv.logit(delta + nu[j]), log = T) +
    dbinom(test_df$Nt[j] - test_df$nt[j], size = M_curr - I, prob = boot::inv.logit(nu[j]), log = T)
}

log_I_posterior_just_rand <- function(I) {
  approxfun(control$I_seq, I_logprior[j, ])(I) +
    approxfun(control$I_seq, I_loglik_just_rand)(I)
}


##################################
# Estimate nu from a fixed point algorithm alternating between 
# updating I and nu 
nu_curr <- nu_approx
converged <- FALSE
print("Running fixed point algorithm")
while (!converged) {
  #####################################################
  # Evaluate the I marginal
  del_fine_by <- .01
  I_fun_list <- I_log_post_eval_list <- list()
  for (I_curr in I_fine) {# I_curr <- I_fine[500]#
    coarse_log_marginal <- log_joint_posterior(delta = del_seq, I = I_curr, del_seq = del_seq, nu = nu_curr)
    if (all(is.na(coarse_log_marginal)) || all(na.omit(coarse_log_marginal) < log_bound)) {
      I_fun_list[[as.character(I_curr)]] <- function(x) -Inf
      I_log_post_eval_list[[as.character(I_curr)]] <- data.frame(I = 0, delta = 3, log_post = -Inf)
    } else {
      coarse_marginal_extreme_delta <- del_seq[range(which(coarse_log_marginal > log_bound))]
      del_seq_fine <- unique(seq(coarse_marginal_extreme_delta[1], coarse_marginal_extreme_delta[2], by = del_fine_by))
      fine_marginal <- log_joint_posterior(delta = del_seq_fine, I = I_curr, del_seq = del_seq_fine, nu = nu_curr)
      log_marg_post_fun <- approxfun(c(0, del_seq_fine, max(del_seq)), c(fine_marginal[1], fine_marginal, fine_marginal[length(fine_marginal)]))
      I_fun_list[[as.character(I_curr)]] <- log_marg_post_fun
      I_log_post_eval_list[[as.character(I_curr)]] <- data.frame(I = I_curr, delta = del_seq_fine, log_post = fine_marginal)
    }
  }
  I_max_log_post <- max(sapply(I_log_post_eval_list, function(x) max(x$log_post)))
  I_marg_post_eval_unnorm <- sapply(I_log_post_eval_list, function(x) sum(exp(x$log_post - I_max_log_post)))
  optimised_I_for_this_nu <- I_fine[which.max(I_marg_post_eval_unnorm)]
  nu_new <- boot::logit((test_df$Nt - test_df$nt) / (test_df$M - optimised_I_for_this_nu))
  if (nu_new - nu_curr < nu_fixed_pt_convergence_tol) {
    converged <- TRUE
  }
  nu_curr <- nu_new
  print(paste("nu =", nu_new))
}
nu_fixed_pt_estimate <- nu_curr
print("Converged")

redo_figures <- FALSE
ltla_df_use <- ltla_df[ltla_df$mid_week == mid_week_use & ltla_df$phe_region == region_curr, ]

bias_store <- list()
for (full_dist_type in c("misspecified", "improved")) {
  nu <- switch(full_dist_type, 
               improved = nu_fixed_pt_estimate, 
               misspecified = nu_approx)
  #####################################################
  # Evaluate the delta marginal
  log_post_I_delta <- matrix(NA, length(del_seq), length(I_seq_coarse), dimnames = list(del_seq, I_seq_coarse))
  fun_list <- log_post_eval_list <- list()
  del_marg <- c()
  for (del_curr in del_seq) {
    coarse_log_marginal <- log_joint_posterior(delta = del_curr, I = I_seq_coarse, del_seq = del_seq, nu = nu)
    log_post_I_delta[as.character(del_curr), ] <- coarse_log_marginal
    if (all(is.na(coarse_log_marginal)) || all(na.omit(coarse_log_marginal) < log_bound)) {
      fun_list[[as.character(del_curr)]] <- function(x) -Inf
      log_post_eval_list[[as.character(del_curr)]] <- data.frame(I = 0, delta = del_curr, log_post = -Inf)
    } else {
      coarse_marginal_extreme_I <- I_seq_coarse[range(which(coarse_log_marginal > log_bound))]
      I_seq_fine <- unique(round(seq(coarse_marginal_extreme_I[1], coarse_marginal_extreme_I[2], by = 100)))
      fine_marginal <- log_joint_posterior(delta = del_curr, I = I_seq_fine, del_seq = del_seq, nu = nu)
      log_marg_post_fun <- approxfun(c(0, I_seq_fine, max(I_seq_coarse)), c(fine_marginal[1], fine_marginal, fine_marginal[length(fine_marginal)]))
      fun_list[[as.character(del_curr)]] <- log_marg_post_fun
      log_post_eval_list[[as.character(del_curr)]] <- data.frame(I = I_seq_fine, delta = del_curr, log_post = fine_marginal)
    }
  }
  max_log_post <- max(sapply(log_post_eval_list, function(x) max(x$log_post)))
  delta_marg_post_eval_unnorm <- sapply(log_post_eval_list, function(x) sum(exp(x$log_post - max_log_post)))
  delta_bin_width <- diff(del_seq[1:2])
  delta_full_marg_post_norm <- delta_marg_post_eval_unnorm / sum(delta_marg_post_eval_unnorm)# / delta_bin_width
  
  plot(del_seq, delta_full_marg_post_norm, ty = "l")
  
  #####################################################
  # Evaluate the bivariate posteriors, joint and cut
  pi_delta_bin_area <- diff(I_fine[1:2]) / M_curr * diff(del_seq[1:2])
  log_post_I_delta_fine <- log_cut_post_I_delta_fine <- matrix(NA, length(del_seq), length(I_fine), dimnames = list(del_seq, I_fine))
  for (del_curr in del_seq) {
    fine_log_marginal <- log_joint_posterior(delta = del_curr, I = I_fine, del_seq = del_seq, nu = nu)
    log_post_I_delta_fine[match(del_curr, del_seq), ] <- fine_log_marginal
  }
  log_post_I_just_rand <- log_I_posterior_just_rand(I_fine)
  log_post_I_just_rand[is.na(log_post_I_just_rand )] <- -Inf
  post_I_just_rand_unnorm <- exp(log_post_I_just_rand - max(log_post_I_just_rand, na.rm = T))
  post_I_just_rand_norm <- post_I_just_rand_unnorm / sum(post_I_just_rand_unnorm, na.rm = T)
  for (I_curr in I_fine) {
    del_log_post_curr <- log_delta_cond_posterior(del_seq, I_curr, del_seq, nu = nu)
    del_post_curr_unnorm <- exp(del_log_post_curr - max(del_log_post_curr, na.rm = T))
    del_post_curr_norm <- del_post_curr_unnorm / sum(del_post_curr_unnorm, na.rm = T)
    log_cut_post_I_delta_fine[, match(I_curr, I_fine)] <- log(del_post_curr_norm) + log(post_I_just_rand_norm[match(I_curr, I_fine)])
  }
  log_post_I_delta_fine[is.na(log_post_I_delta_fine)] <- -Inf
  log_cut_post_I_delta_fine[is.na(log_cut_post_I_delta_fine)] <- -Inf
  joint_2d_log_posterior_unnorm <- log_post_I_delta_fine - max(log_post_I_delta_fine)
  joint_2d_posterior_unnorm <- exp(joint_2d_log_posterior_unnorm)
  joint_2d_posterior_norm <- joint_2d_posterior_unnorm / sum(joint_2d_posterior_unnorm, na.rm = T)# / pi_delta_bin_area
  joint_2d_log_posterior <- log(joint_2d_posterior_norm)
  
  joint_2d_cut_posterior_unnorm <- exp(log_cut_post_I_delta_fine - max(log_cut_post_I_delta_fine))
  joint_2d_cut_posterior_norm <- joint_2d_cut_posterior_unnorm / sum(joint_2d_cut_posterior_unnorm)# / pi_delta_bin_area
  joint_2d_log_cut_posterior <- log(joint_2d_cut_posterior_norm)
  delta_cut_marg_post_norm <- rowSums(joint_2d_cut_posterior_norm)
  pi_cut_marg_post_norm <- colSums(joint_2d_cut_posterior_norm)
  
  
  #####################################################
  # Evaluate the I marginal
  del_fine_by <- .01
  I_fun_list <- I_log_post_eval_list <- list()
  for (I_curr in I_fine) {# I_curr <- I_fine[500]#
    coarse_log_marginal <- log_joint_posterior(delta = del_seq, I = I_curr, del_seq = del_seq, nu = nu)
    if (all(is.na(coarse_log_marginal)) || all(na.omit(coarse_log_marginal) < log_bound)) {
      I_fun_list[[as.character(I_curr)]] <- function(x) -Inf
      I_log_post_eval_list[[as.character(I_curr)]] <- data.frame(I = 0, delta = del_curr, log_post = -Inf)
    } else {
      coarse_marginal_extreme_delta <- del_seq[range(which(coarse_log_marginal > log_bound))]
      del_seq_fine <- unique(seq(coarse_marginal_extreme_delta[1], coarse_marginal_extreme_delta[2], by = del_fine_by))
      fine_marginal <- log_joint_posterior(delta = del_seq_fine, I = I_curr, del_seq = del_seq_fine, nu = nu)
      log_marg_post_fun <- approxfun(c(0, del_seq_fine, max(del_seq)), c(fine_marginal[1], fine_marginal, fine_marginal[length(fine_marginal)]))
      I_fun_list[[as.character(I_curr)]] <- log_marg_post_fun
      I_log_post_eval_list[[as.character(I_curr)]] <- data.frame(I = I_curr, delta = del_seq_fine, log_post = fine_marginal)
    }
  }
  I_max_log_post <- max(sapply(I_log_post_eval_list, function(x) max(x$log_post)))
  I_marg_post_eval_unnorm <- sapply(I_log_post_eval_list, function(x) sum(exp(x$log_post - I_max_log_post)))
  pi_bin_width <- diff(pi_seq_coarse[1:2])
  pi_full_marg_post_norm <- I_marg_post_eval_unnorm / sum(I_marg_post_eval_unnorm)# / pi_bin_width
  
  
  I_quant <- prevdebiasr:::randomised_testing_prevalence(test_df, control, imperfect)
  delta_regional <- prevdebiasr:::delta_regional_posterior(test_df, I_quant, control, imperfect)
  delta_df_cut_for_inference <- data.frame(delta_prior_mean = delta_regional$delta_post_mean,
                                           delta_prior_sd = delta_regional$delta_post_sd)
  
  delta_post_moment1 <- sum(delta_full_marg_post_norm * del_seq, na.rm = TRUE)
  delta_post_moment2 <- sum(delta_full_marg_post_norm * del_seq^2, na.rm = TRUE)
  delta_post_mean <- delta_post_moment1
  delta_post_sd <- sqrt(delta_post_moment2 - (delta_post_moment1)^2)
  delta_df_full_for_inference <- data.frame(delta_prior_mean = delta_post_mean,
                                              delta_prior_sd = delta_post_sd)
  
  
  #####################################################
  # Evaluate the LTLA-level prevalence for joint and cut
  ltla_prev_list <- list()
  for(plot_type in c("full", "cut")) {
    if(plot_type == "full") {
      delta_df <- cbind(test_df[j, c("phe_region", "mid_week")], delta_df_full_for_inference)
    }
    if(plot_type == "cut") {
      delta_df <- cbind(test_df[j, c("phe_region", "mid_week")], delta_df_cut_for_inference)
    }
    ltla_list <- ltla_df_use %>%
      left_join(delta_df, by = c("phe_region", "mid_week")) %>%
      group_by(ltla) %>%
      group_split()
    ltla_names <- sapply(ltla_list, function(x) x$ltla[1])
    ltla_prevalence_curr <- parLapply(clust, ltla_list, local_prevalence,
                                 control, imperfect, type)
    names(ltla_prevalence_curr) <- ltla_names
    postmat <- sapply(ltla_prevalence_curr, function(x) x$norm_post)
    postmat_cum <- apply(postmat, 2, cumsum)
    ci_lower <- control$I_seq[apply(postmat_cum, 2, function(x) findInterval(.025, x))]
    ci_upper <- control$I_seq[apply(postmat_cum, 2, function(x) findInterval(.975, x))]
    ci_mid <- control$I_seq[apply(postmat_cum, 2, function(x) findInterval(.5, x))]
    mom1 <- colSums(postmat * control$I_seq)
    mom2 <- colSums(postmat * control$I_seq^2)
    sdv <- sqrt(mom2 - mom1^2)
    ltla_df_use[match(ltla_names, ltla_df_use$ltla), paste0(plot_type, "_prevnum_post_mn")] <- ci_mid
    ltla_df_use[match(ltla_names, ltla_df_use$ltla), paste0(plot_type, "_prevnum_post_sd")] <- sdv
    ltla_df_use[match(ltla_names, ltla_df_use$ltla), paste0(plot_type, "_prevnum_post_lower")] <- ci_lower
    ltla_df_use[match(ltla_names, ltla_df_use$ltla), paste0(plot_type, "_prevnum_post_upper")] <- ci_upper
    ltla_df_use[, paste0(plot_type, "_prevprop_post_mn")] <- ltla_df_use[, paste0(plot_type, "_prevnum_post_mn")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_sd")] <- ltla_df_use[, paste0(plot_type, "_prevnum_post_sd")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_lower")] <- ltla_df_use[, paste0(plot_type, "_prevnum_post_lower")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_upper")] <- ltla_df_use[, paste0(plot_type, "_prevnum_post_upper")] / ltla_df_use$M
    ltla_prev_list[[plot_type]] <- ltla_prevalence_curr
  }
  
  
  
  #####################################################
  # Prepare gold standard react LTLA round-aggregated estimates
  react_date_df <- readr::read_csv(file.path(data_dir_to_use, "react_dates.csv"))
  react_ltla_df <- readr::read_csv(file.path(data_dir_to_use, "react_ltla.csv"))
  this_round <- "8"
  comp_1 <- react_ltla_df %>%
    ungroup() %>%
    filter(round == this_round) %>%
    select(ltla, l, m, u)
  comp_1 <- comp_1[match(ltla_df_use$ltla, comp_1$ltla), ]
  
  sum(exp(joint_2d_log_cut_posterior))
  
  
  
  if(nu == nu_approx) {
    joint_2d_log_posterior_misspecified <- joint_2d_log_posterior
    delta_full_marg_post_norm_misspecified <- delta_full_marg_post_norm
    pi_full_marg_post_norm_misspecified <- pi_full_marg_post_norm
    delta_df_full_for_inference_misspecified <- data.frame(delta_prior_mean = delta_post_mean,
                                              delta_prior_sd = delta_post_sd)
    plot_type <- "full_misspecified"
    ltla_df_use[, paste0(plot_type, "_prevprop_post_mn")] <- ltla_df_use[, paste0("full_prevnum_post_mn")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_sd")] <- ltla_df_use[, paste0("full_prevnum_post_sd")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_lower")] <- ltla_df_use[, paste0("full_prevnum_post_lower")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_upper")] <- ltla_df_use[, paste0("full_prevnum_post_upper")] / ltla_df_use$M
    
  }
  if(nu == nu_fixed_pt_estimate) {
    joint_2d_log_posterior_improved <- joint_2d_log_posterior
    delta_full_marg_post_norm_improved <- delta_full_marg_post_norm
    pi_full_marg_post_norm_improved <- pi_full_marg_post_norm
    delta_df_full_for_inference_improved <- data.frame(delta_prior_mean = delta_post_mean,
                                                           delta_prior_sd = delta_post_sd)
    plot_type <- "full_improved"
    ltla_df_use[, paste0(plot_type, "_prevprop_post_mn")] <- ltla_df_use[, paste0("full_prevnum_post_mn")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_sd")] <- ltla_df_use[, paste0("full_prevnum_post_sd")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_lower")] <- ltla_df_use[, paste0("full_prevnum_post_lower")] / ltla_df_use$M
    ltla_df_use[, paste0(plot_type, "_prevprop_post_upper")] <- ltla_df_use[, paste0("full_prevnum_post_upper")] / ltla_df_use$M
  }
  
  delta_df_full_for_inference <- data.frame(delta_prior_mean = delta_post_mean,
                                            delta_prior_sd = delta_post_sd)
  
  if (redo_figures) {
    #####################################################
    # Plot it!
    grey_pallette <- grey(seq(0, 1, length.out = 1000))
    graphics.off()
    # plot_dir <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/figures"
    plot_file <- paste0(plot_dir_curr, "/cut_vs_full", ifelse(nu == nu_approx, "", "_nu_fixed_pt"), ".jpeg")
    jpeg(plot_file, 9, 9, res = 750, units = "in")
    par(mar = c(3, 3, 5, 5), oma = c(1, 1, 1, 1), mfrow = c(2, 2))
    cexax <- 1
    # plot_bd <- -20
    c("full_misspecified", "cut", "full_fixed")
    for(plot_type in c("full", "cut")) {
      if(plot_type == "full") {
        zpl <- joint_2d_log_posterior
        x_marg <- delta_full_marg_post_norm
        y_marg <- pi_full_marg_post_norm
      }
      if(plot_type == "cut") {
        zpl <- joint_2d_log_cut_posterior
        x_marg <- delta_cut_marg_post_norm
        y_marg <- pi_cut_marg_post_norm
      }
      cut_diff_from_max_log_post <- 20
      # cut_at_for_plot <- floor((max(c(joint_2d_log_posterior, joint_2d_log_cut_posterior)) - cut_diff_from_max_log_post) / 10) * 10
      max_raw_log_post <- max(c(joint_2d_log_posterior, joint_2d_log_cut_posterior))
      cut_at_for_plot <- max_raw_log_post - cut_diff_from_max_log_post
      raw_range_for_plot <- c(cut_at_for_plot, max_raw_log_post)
      # bin_area_corrected_range_for_legend <- raw_range_for_plot - log(pi_delta_bin_area)
      # zpl[zpl < cut_at_for_plot] <- cut_at_for_plot
      zpl_scaled <- (zpl - raw_range_for_plot[1]) / diff(raw_range_for_plot)
      zpl_scaled[zpl_scaled < 0] <- 0
      # TODO: scale bar for log posteriors, tricky thing being getting the "< x" lower bound to be a round number, 
      # when accounting for the bin area to give the true log posterior density (as opposed to normalised by bin)
      
      image(x = 1:n_del, y = 1:n_pi_fine_plot, z = zpl_scaled, xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
            col = grey_pallette, bty = "n", zlim = c(0, 1))
      mtext(side = 2, line = 2.5, text = expression(paste("% Prevalence  ", pi)), cex = cexax)
      mtext(side = 1, line = 2.5, text = expression(paste("Bias parameter  ", delta)), cex = cexax)
      del_labs <- pretty(del_seq)
      del_labs <- del_labs[del_labs > min(del_seq) & del_labs < max(del_seq)]
      del_ats <- findInterval(del_labs, del_seq)
      axis(side = 1, at = del_ats, labels = del_labs, las = 0)
      pi_seq_plot <- I_fine / M_curr * 100
      pi_labs <- pretty(pi_seq_plot)
      pi_labs <- pi_labs[pi_labs / 100 * M_curr >= min(I_fine) & pi_labs / 100 * M_curr <= max(I_fine)]
      pi_ats <- findInterval(pi_labs, pi_seq_plot) + 1
      axis(side = 2, at = pi_ats, labels = pi_labs, las = 2)
      par(xpd = NA)
      yran <- c(1, n_pi_fine_plot)
      xran <- c(1, n_del)
      rel_hei <- .1
      new_y <- yran[2] + x_marg * .1 * diff(yran) / max(x_marg)
      new_x <- xran[2] + y_marg * .1 * diff(xran) / max(y_marg)
      lines(1:n_del, new_y)
      lines(new_x, 1:n_pi_fine_plot)
      del_mean_curr <- switch (plot_type,
                               full = delta_df_full_for_inference$delta_prior_mean,
                               cut = delta_df_cut_for_inference$delta_prior_mean
      )
      I_mean_curr <- sum(y_marg * I_fine)
      x_loc <- findInterval(del_mean_curr, del_seq)
      y_loc <- findInterval(I_mean_curr, I_fine)
      lines(rep(x_loc, 2), range(new_y), lty = 3, col = 1)
      lines(range(new_x), rep(y_loc, 2), lty = 3, col = 1)
      text(x_loc, range(new_y)[2] + diff(range(new_y)) * .4, formatC(del_mean_curr, format = "f", digits = 1), adj = c(0.5, 1), cex = .8)
      text(range(new_x)[2] + diff(range(new_x)) * .75, y_loc, formatC(I_mean_curr / M_curr * 100, format = "f", digits = 1), adj = c(1, .5), cex = .9)
      par(xpd = F)
      mtext_curr <- switch (plot_type,
                            full = "Full posterior",
                            cut = "Cut posterior"
      )
      mtext(side = 3, line = 3.5, text = mtext_curr, cex = 1.2)
      mtext(side = 3, line = 1, at = 0, text = switch (plot_type, full = "(a)", cut = "(b)"), cex = 1)
    }
    for (plot_type in c("full", "cut")) {
      comp_2 <- data.frame(ltla = ltla_df_use$ltla, 
                           l = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_lower")]) * 100,
                           m = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_mn")]) * 100, 
                           u = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_upper")]) * 100)
      plot(comp_1$m, comp_2$m, xlim = c(0, 5), ylim = c(0, 5), ty = "n", xlab = "", ylab = "", las = 1, xaxs = "i", yaxs = "i")
      mtext(side = 1, line = 2.5, text = "% Prevalence (REACT)", cex = cexax)
      mtext(side = 2, line = 2.5, text = "% Prevalence (debiased Pillar 1+2)", cex = cexax)
      bias_store_name <- ifelse(plot_type == "full", 
                                ifelse(nu == nu_approx, "full_misspecified", "full_improved"),
                                "cut")
      bias_store[[paste0(bias_store_name, "_mean")]] <- bias_mean <- mean(comp_2$m - comp_1$m)
      bias_store[[paste0(bias_store_name, "_se")]] <- bias_se <- sd(comp_2$m - comp_1$m) / sqrt(nrow(comp_1))
      # mtext(side = 3, line = 0.25, text = paste0("Bias = ", formatC(x = bias_mean, format = "f", digits = 2),
      #                                            "% (SE = ", formatC(x = bias_se, format = "f", digits = 2), "%)"), cex = cexax)
      abline(0, 1)
      for (k in 1:nrow(comp_1)) {
        lines(x = rep(comp_1$m[k], 2), y = unlist(comp_2[k, c("l", "u")]))
        lines(x = unlist(comp_1[k, c("l", "u")]), y = rep(comp_2$m[k], 2), col = "grey")
      }
      points(comp_1$m, comp_2$m, pch = 19, cex = .7)
      mtext(side = 3, line = 1, at = 0, text = switch (plot_type, full = "(c)", cut = "(d)"), cex = 1)
    }
    dev.off()
  } else {
    for (plot_type in c("full", "cut")) {
      comp_2 <- data.frame(ltla = ltla_df_use$ltla, 
                           l = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_lower")]) * 100,
                           m = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_mn")]) * 100, 
                           u = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_upper")]) * 100)
      bias_store_name <- ifelse(plot_type == "full", 
                                ifelse(nu == nu_approx, "full_misspecified", "full_improved"),
                                "cut")
      bias_store[[paste0(bias_store_name, "_mean")]] <- bias_mean <- mean(comp_2$m - comp_1$m)
      bias_store[[paste0(bias_store_name, "_se")]] <- bias_se <- sd(comp_2$m - comp_1$m) / sqrt(nrow(comp_1))
    }  
  }
}  

# Output means and 95% CIs for marg delta
delta_mean_cut <- sum(delta_cut_marg_post_norm * del_seq)
delta_sd_cut <- delta_df_cut_for_inference$delta_prior_sd
delta_lower_cut <- del_seq[findInterval(0.025, cumsum(delta_cut_marg_post_norm)) - 1]
delta_upper_cut <- del_seq[findInterval(0.975, cumsum(delta_cut_marg_post_norm)) - 1]

delta_mean_full_improved <- sum(delta_full_marg_post_norm_improved * del_seq)
delta_sd_full_improved <- delta_df_full_for_inference_improved$delta_prior_sd
delta_lower_full_improved <- del_seq[findInterval(0.025, cumsum(delta_full_marg_post_norm_improved)) - 1]
delta_upper_full_improved <- del_seq[findInterval(0.975, cumsum(delta_full_marg_post_norm_improved)) - 1]

delta_mean_full_misspecified <- sum(delta_full_marg_post_norm_misspecified * del_seq)
delta_sd_full_misspecified <- delta_df_full_for_inference_misspecified$delta_prior_sd
delta_lower_full_misspecified <- del_seq[findInterval(0.025, cumsum(delta_full_marg_post_norm_misspecified)) - 1]
delta_upper_full_misspecified <- del_seq[findInterval(0.975, cumsum(delta_full_marg_post_norm_misspecified)) - 1]

delta_full_marg_post_norm_improved


pi_mean_cut <- sum(pi_cut_marg_post_norm * pi_seq_coarse) * 100
pi_lower_cut <- pi_seq_coarse[findInterval(0.025, cumsum(pi_cut_marg_post_norm)) - 1] * 100
pi_upper_cut <- pi_seq_coarse[findInterval(0.975, cumsum(pi_cut_marg_post_norm)) - 1] * 100

pi_mean_full_improved <- sum(pi_full_marg_post_norm_improved * pi_seq_coarse) * 100
pi_lower_full_improved <- pi_seq_coarse[findInterval(0.025, cumsum(pi_full_marg_post_norm_improved)) - 1] * 100
pi_upper_full_improved <- pi_seq_coarse[findInterval(0.975, cumsum(pi_full_marg_post_norm_improved)) - 1] * 100

pi_mean_full_misspecified <- sum(pi_full_marg_post_norm_misspecified * pi_seq_coarse) * 100
pi_lower_full_misspecified <- pi_seq_coarse[findInterval(0.025, cumsum(pi_full_marg_post_norm_misspecified)) - 1] * 100
pi_upper_full_misspecified <- pi_seq_coarse[findInterval(0.975, cumsum(pi_full_marg_post_norm_misspecified)) - 1] * 100


bias_mean_full_improved <- bias_store$full_improved_mean
bias_se_full_improved <- bias_store$full_improved_se

bias_mean_full_misspecified <- bias_store$full_misspecified_mean
bias_se_full_misspecified <- bias_store$full_misspecified_se

bias_mean_cut <- bias_store$cut_mean
bias_se_cut <- bias_store$cut_se

# dir_text_numbers <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/text_numbers/cut_vs_full_comp"
save_num <- c("delta_mean_cut", "delta_sd_cut", "delta_lower_cut", "delta_upper_cut", 
              "delta_mean_full_improved", "delta_sd_full_improved", "delta_lower_full_improved", "delta_upper_full_improved", 
              "delta_mean_full_misspecified", "delta_sd_full_misspecified", "delta_lower_full_misspecified", "delta_upper_full_misspecified", 
              "pi_mean_cut", "pi_lower_cut", "pi_upper_cut", 
              "pi_mean_full_misspecified", "pi_lower_full_misspecified", "pi_upper_full_misspecified",
              "pi_mean_full_improved", "pi_lower_full_improved", "pi_upper_full_improved")


for(numc in save_num)
  write.table(formatC(eval(as.name(numc)), format = "f", digits = 1), file = paste(dir_text_numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)

save_num <- c("nu_fixed_pt_estimate", "nu_approx")
for(numc in save_num)
  write.table(formatC(eval(as.name(numc)), format = "f", digits = 3), file = paste(dir_text_numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


save_num <- c("bias_mean_full_improved", "bias_se_full_improved", 
              "bias_mean_full_misspecified", "bias_se_full_misspecified", 
              "bias_mean_cut", "bias_se_cut")
for(numc in save_num)
  write.table(formatC(eval(as.name(numc)), format = "f", digits = 2), file = paste(dir_text_numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)







#####################################################
# Plot it!
grey_pallette <- grey(seq(0, 1, length.out = 1000))
graphics.off()
# plot_dir <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/figures"
plot_file <- paste0(plot_dir_curr, "/cut_vs_full_3-way_with_improved_nu_estimate.jpeg")
jpeg(plot_file, 12, 8, res = 750, units = "in")
par(mar = c(3, 3, 5, 5), oma = c(1, 1, 3, 1), mfrow = c(2, 3))
cexax <- 1
# plot_bd <- -20
for(plot_type in c("full_misspecified", "cut", "full_improved")) {
  if(plot_type == "full_misspecified") {
    zpl <- joint_2d_log_posterior_misspecified
    x_marg <- delta_full_marg_post_norm_misspecified
    y_marg <- pi_full_marg_post_norm_misspecified
  }
  if(plot_type == "cut") {
    zpl <- joint_2d_log_cut_posterior
    x_marg <- delta_cut_marg_post_norm
    y_marg <- pi_cut_marg_post_norm
  }
  if(plot_type == "full_improved") {
    zpl <- joint_2d_log_posterior_improved
    x_marg <- delta_full_marg_post_norm_improved
    y_marg <- pi_full_marg_post_norm_improved
  }
  cut_diff_from_max_log_post <- 20
  # cut_at_for_plot <- floor((max(c(joint_2d_log_posterior, joint_2d_log_cut_posterior)) - cut_diff_from_max_log_post) / 10) * 10
  max_raw_log_post <- max(c(joint_2d_log_posterior_misspecified, joint_2d_log_cut_posterior, joint_2d_log_posterior_improved))
  cut_at_for_plot <- max_raw_log_post - cut_diff_from_max_log_post
  raw_range_for_plot <- c(cut_at_for_plot, max_raw_log_post)
  # bin_area_corrected_range_for_legend <- raw_range_for_plot - log(pi_delta_bin_area)
  # zpl[zpl < cut_at_for_plot] <- cut_at_for_plot
  zpl_scaled <- (zpl - raw_range_for_plot[1]) / diff(raw_range_for_plot)
  zpl_scaled[zpl_scaled < 0] <- 0
  # TODO: scale bar for log posteriors, tricky thing being getting the "< x" lower bound to be a round number, 
  # when accounting for the bin area to give the true log posterior density (as opposed to normalised by bin)
  
  image(x = 1:n_del, y = 1:n_pi_fine_plot, z = zpl_scaled, xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
        col = grey_pallette, bty = "n", zlim = c(0, 1))
  mtext(side = 2, line = 2.5, text = expression(paste("% Prevalence  ", pi)), cex = cexax)
  mtext(side = 1, line = 2.5, text = expression(paste("Bias parameter  ", delta)), cex = cexax)
  del_labs <- pretty(del_seq)
  del_labs <- del_labs[del_labs > min(del_seq) & del_labs < max(del_seq)]
  del_ats <- findInterval(del_labs, del_seq)
  axis(side = 1, at = del_ats, labels = del_labs, las = 0)
  pi_seq_plot <- I_fine / M_curr * 100
  pi_labs <- pretty(pi_seq_plot)
  pi_labs <- pi_labs[pi_labs / 100 * M_curr >= min(I_fine) & pi_labs / 100 * M_curr <= max(I_fine)]
  pi_ats <- findInterval(pi_labs, pi_seq_plot) + 1
  axis(side = 2, at = pi_ats, labels = pi_labs, las = 2)
  par(xpd = NA)
  yran <- c(1, n_pi_fine_plot)
  xran <- c(1, n_del)
  rel_hei <- .1
  new_y <- yran[2] + x_marg * .1 * diff(yran) / max(x_marg)
  new_x <- xran[2] + y_marg * .1 * diff(xran) / max(y_marg)
  lines(1:n_del, new_y)
  lines(new_x, 1:n_pi_fine_plot)
  del_mean_curr <- switch (plot_type,
                           full_misspecified = delta_df_full_for_inference_misspecified$delta_prior_mean,
                           cut = delta_df_cut_for_inference$delta_prior_mean,
                           full_improved = delta_df_full_for_inference_improved$delta_prior_mean,
  )
  I_mean_curr <- sum(y_marg * I_fine)
  x_loc <- findInterval(del_mean_curr, del_seq)
  y_loc <- findInterval(I_mean_curr, I_fine)
  lines(rep(x_loc, 2), range(new_y), lty = 3, col = 1)
  lines(range(new_x), rep(y_loc, 2), lty = 3, col = 1)
  text(x_loc, range(new_y)[2] + diff(range(new_y)) * .4, formatC(del_mean_curr, format = "f", digits = 1), adj = c(0.5, 1), cex = .95)
  text(range(new_x)[2] + diff(range(new_x)) * .75, y_loc, formatC(I_mean_curr / M_curr * 100, format = "f", digits = 1), adj = c(1, .5), cex = .95)
  par(xpd = F)
  mtext_curr_outer <- switch (plot_type,
                              full_misspecified = "Full posterior",
                              cut = "Cut posterior",
                              full_improved = "Full posterior"
  )
  mtext(side = 3, line = 5.5, text = mtext_curr_outer, cex = 1.2)
  # mtext_curr_inner <- switch (plot_type,
  #                             full_misspecified = "Misspecified model",
  #                             cut = "Misspecified model",
  #                             full_improved = "Improved model"
  # )
  # mtext(side = 3, line = 3.5, text = mtext_curr_inner, cex = .85)
  cex_inner <- 1
  line_inner <- 3.5
  if (plot_type == "full_misspecified") {
    mtext(side = 3, line = line_inner, text = expression(paste("Misspecified model,", nu == hat(nu))), cex = cex_inner)
  }
  if (plot_type == "cut") {
    mtext(side = 3, line = line_inner, text = expression(paste("Misspecified model,", nu == hat(nu))), cex = cex_inner)
  }
  if (plot_type == "full_improved") {
    mtext(side = 3, line = line_inner, text = expression(paste("Improved model,", nu == hat(nu)[MLE])), cex = cex_inner)
  }
  mtext(side = 3, line = 1, at = 0, text = switch (plot_type, full_misspecified = "(a)", cut = "(b)", full_improved = "(c)"), cex = 1)
}
bias_store <- list()
for (plot_type in c("full_misspecified", "cut", "full_improved")) {
  comp_2 <- data.frame(ltla = ltla_df_use$ltla, 
                       l = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_lower")]) * 100,
                       m = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_mn")]) * 100, 
                       u = unlist(ltla_df_use[, paste0(plot_type, "_prevprop_post_upper")]) * 100)
  plot(comp_1$m, comp_2$m, xlim = c(0, 5), ylim = c(0, 5), ty = "n", xlab = "", ylab = "", las = 1, xaxs = "i", yaxs = "i")
  mtext(side = 1, line = 2.5, text = "% Prevalence (REACT)", cex = cexax)
  mtext(side = 2, line = 2.5, text = "% Prevalence (debiased Pillar 1+2)", cex = cexax)
  bias_store[[paste0(plot_type, "_mean")]] <- bias_mean <- mean(comp_2$m - comp_1$m)
  bias_store[[paste0(plot_type, "_se")]] <- bias_se <- sd(comp_2$m - comp_1$m) / sqrt(nrow(comp_1))
  # mtext(side = 3, line = 0.25, text = paste0("Bias = ", formatC(x = bias_mean, format = "f", digits = 2),
  #                                            "% (SE = ", formatC(x = bias_se, format = "f", digits = 2), "%)"), cex = cexax)
  abline(0, 1)
  for (k in 1:nrow(comp_1)) {
    lines(x = rep(comp_1$m[k], 2), y = unlist(comp_2[k, c("l", "u")]))
    lines(x = unlist(comp_1[k, c("l", "u")]), y = rep(comp_2$m[k], 2), col = "grey")
  }
  points(comp_1$m, comp_2$m, pch = 19, cex = .7)
  mtext(side = 3, line = 1, at = 0, text = switch (plot_type, full_misspecified = "(d)", cut = "(e)", full_improved = "(f)"), cex = 1)
}
dev.off()


























