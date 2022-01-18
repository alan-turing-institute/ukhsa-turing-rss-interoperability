############################## 
### SIR control parameters
############################## 
control_SIR <- list()

# Temporal autocorrelation prior on R number (determines the transition matrix for HMM on beta)
control_SIR$R_AR_sd <- 0.2

# Prior on proportion of immune population at t = 0
control_SIR$prior_immune_t0 <- list(mean = 0.06,
                                    sd = 0.01,
                                    upper_trunc = 0.1,
                                    lower_trunc = 0)

# MCMC parameters 
control_SIR$n_iters <- 200
control_SIR$burn_in <- 100

# # MCMC parameters for checking
# control_SIR$n_iters <- 600
# control_SIR$burn_in <- 300


# R_grid is the grid for R number (this determines the discrete state space of the HMM on beta)
# Note that R = beta / gamma, so the grid on R determines the grid on beta (as gamma is prespecified)
control_SIR$R_grid <- seq(.05, 4, by = .05)

# Parameters for epidemic model
control_SIR$generation_time <- 1 # Epidemic generation time parameter - recovery rate \gamma = 1 / generation_time
control_SIR$epi_gamma <- 1 / control_SIR$generation_time

############################################
### Functions to fit SIR epidemic model
############################################ 

# Main wrapper function

# d = d_ltla

sample_I_R_SIR_from_pre_calc_lik <- function(d, 
                                             I_log_lik, 
                                             trans_mats,
                                             control_debias,
                                             control_SIR, 
                                             init_meth = c("biased_data", 
                                                           "random_data")[1]) {
  if (all(c("ltla", "mid_week") %in% colnames(d))) {
    d$namdat <- paste0(d$ltla, "_", d$mid_week)
  } else {
    d$namdat <- 1:nrow(d)
  }
  d$I <- control_debias$I_seq[apply(I_log_lik, 1, function(v){ if(all(is.na(v))) {NA} else {which.max(v)}})]
  
  d$I <- c(mean(d$I[1:2]),
           stats::filter(d$I, rep(1 / 3, 3))[-c(1, nrow(d))],
           mean(d$I[-(1:(nrow(d) - 2))]))
  d$I <- control_debias$I_seq[findInterval(d$I, control_debias$I_seq)]
  d$dI <- c(diff(d$I), 0)
  R_rough_estimate <- (d$dI + d$I * control_SIR$epi_gamma) / pmax(d$I, 1) / control_SIR$epi_gamma
  R_trans <- sapply(trans_mats, function(x) {
    out <- c()
    for (i in 1:(nrow(d) - 1)) {
      out <- c(out, x[as.character(d$I[i]), as.character(d$I[i+1])])
    }
    out
  })
  d$R <- c(control_SIR$R_grid[apply(R_trans, 1, which.max)], 1)
  saml <- rep(times = 3, x = list(matrix(NA, nrow(d), control_SIR$n_iters, 
                                         dimnames = list(d$namdat, 1:control_SIR$n_iters))))
  names(saml) <- c("I", "R", "Rplus")
  itc <- 0
  nI <- length(control_debias$I_seq)
  parl.zero <- parl.inf <- list()
  for (j in 1:nrow(d)) {
    parl.zero[[j]] <- matrix(0, nI, nI)
    parl.inf[[j]] <- matrix(-Inf, nI, nI)
  }
  while (itc < control_SIR$n_iters) {
    itc <- itc + 1
    print(itc)
    try({
      d$I <- sample_I_SIR_from_pre_calc_lik(d, control_debias$I_seq, trans_mats, 
                                            log_lik = I_log_lik,
                                            parl.zero = parl.zero, 
                                            parl.inf = parl.inf)
    })
    d$R <- sample_R_SIR_from_pre_calc_lik(d, control_debias$I_seq, 
                                          control_SIR$R_grid, 
                                          control_SIR$R_AR_sd, trans_mats)
    d$Rplus <- sample_Rplus(d, control_SIR$prior_immune_t0, 
                            control_SIR$epi_gamma)
    saml$I[, itc] <- d$I 
    saml$R[, itc] <- d$R
    saml$Rplus[, itc] <- d$Rplus
  }
  return(saml)
}


# Sample I (infectious) conditional on rest of parameters

sample_I_SIR_from_pre_calc_lik <- function(d, I_seq, trans_mats, 
                                           log_lik, parl.zero = NULL, 
                                           parl.inf = NULL) {
  lik <- exp(log_lik - apply(log_lik, 1, max))
  d$I_index <- NA
  parl <- parl.zero
  parll <- parl.inf
  nI <- length(I_seq)
  if (is.null(parl.zero) | is.null(parl.inf)) {
    parl.zero <- parl.inf <- list()
    for (j in 1:nrow(d)) {
      parl.zero[[j]] <- matrix(0, nI, nI)
      parl.inf[[j]] <- matrix(-Inf, nI, nI)
    }
  }
  # for(j in 1:23){
  for(j in 1:nrow(d)){
      if(j == 1){
      for (k in 1:nI) {
        parl[[j]][k, ] <- lik[j, ]
      }
      parl[[j]] <- parl[[j]] / sum(parl[[j]])
    }
    if (j > 1) {
      prev_fwd_probs <- colSums(parl[[j - 1]])
      updatec <- which(prev_fwd_probs != 0)
      
      parll[[j]][updatec, ] <- sweep(log(prev_fwd_probs[updatec]) + log(trans_mats[[as.character(d$R[j])]][updatec, , drop = F]), 2, log_lik[j, ], '+')
      parl[[j]][updatec, ] <- exp(parll[[j]][updatec, , drop = F] - max(parll[[j]][updatec, , drop = F]))
      parl[[j]][updatec, ] <- parl[[j]][updatec, , drop = F] / sum(parl[[j]][updatec, , drop = F])
      if(length(updatec) == 0) {
        stop()
      }
    }
  }
  for (j in nrow(d):1) {
    if (j == nrow(d)) {
      multinom_sam_probs <- colSums(parl[[j]])
    }
    if (j < nrow(d)) {
      multinom_sam_probs <- parl[[j + 1]][, d$I_index[j + 1]]
    }
    multinom_sam_probs <- multinom_sam_probs / sum(multinom_sam_probs)
    d$I_index[j] <- which(rmultinom(1, 1, multinom_sam_probs) == 1)
  }
  I_out <- I_seq[d$I_index] 
  return(I_out)
}


# I_seq <- control_debias$I_seq
# R_grid <- control_SIR$R_grid 
# R_AR_sd <- control_SIR$R_AR_sd
# trll <- trans_mats
# 
# par(mfrow = c(2, 1))
# plot(d$I)
# plot(d$R)


# Sample effective R conditional on rest of parameters
sample_R_SIR_from_pre_calc_lik <- function(d, I_seq, R_grid, R_AR_sd, trll) {
  
  d$I_index <- match(d$I, I_seq)
  nR <- length(R_grid)
  Rarr <- array(NA, dim = c(nrow(d), nR, nR), 
                dimnames = list(d$namdat, R_grid, R_grid))
  Rtrans <- matrix(NA, nR, nR, 
                   dimnames = list(R_grid, R_grid))
  Rlikmat <- matrix(-Inf, nrow(d), nR, 
                    dimnames = list(d$namdat, R_grid))
  for(Rc in R_grid){
    Rtrans[as.character(Rc), ] <- dnorm(R_grid, Rc, sd = R_AR_sd)
    Rtrans[as.character(Rc), ] <- Rtrans[as.character(Rc), ] / sum(Rtrans[as.character(Rc), ])
  }
  for (j in 1:nrow(d)) {
    if (j == 1) {
      for(k in 1:nR) {
        Rarr[j, k, ] <- 1
      }
    }
    if (j > 1) {
      for (Rc in R_grid) {
        Rlikmat[j, as.character(Rc)] <- trll[[as.character(Rc)]][d$I_index[j - 1], d$I_index[j]]
      }
      Rarr[j, , ] <- sweep(colSums(Rarr[j - 1, , ]) * Rtrans, 2, Rlikmat[j, ], '*')
    }
    if(all(Rarr[j, , ] == 0))
      stop()
    Rarr[j, , ] <- Rarr[j, , ] / sum(Rarr[j, , ])
  }
  for(j in nrow(d):1){
    if(j == nrow(d)){ #  || d$name[j] != d$name[j + 1]
      multpc <- colSums(Rarr[j, , ])
    }
    if(j < nrow(d)){ #  && d$name[j] == d$name[j + 1]
      multpc <- Rarr[j + 1, , as.character(d$R[j + 1])]
    }
    multpc <- multpc / sum(multpc)
    d$R[j] <- R_grid[which(rmultinom(1, 1, multpc) == 1)]
  }
  return(d$R)
}



# Function for sampling immune population Rplus

sample_Rplus <- function(d, prior_immune_t0, epi_gamma, 
                         n.quantile.approx = 25) {
  
  if (is.null(d$V)) d$V <- 0
  
  M.curr.region <- d$M[1]
  d$delRplus <- NA
  quse <- 1:n.quantile.approx / (n.quantile.approx + 1)
  for (j in 1:nrow(d)) {
    if (j == 1) {
      delRplus.dist <- round(truncnorm::qtruncnorm(quse, prior_immune_t0$lower_trunc,
                                                   prior_immune_t0$upper_trunc,
                                                   prior_immune_t0$mean,
                                                   prior_immune_t0$sd) * M.curr.region)
    } else {
      delR.dist <- qbinom(p = quse, size = d$I[j - 1], 
                          prob = stats::pexp(q = 1, rate = epi_gamma))
      delVtilde.dist <- qhyper(p = quse, M.curr.region - d$V[j - 1], 
                               M.curr.region - d$Rplus[j - 1] - d$I[j - 1], 
                               d$V[j] - d$V[j - 1])
      delRplus.dist <- c(outer(delR.dist, delVtilde.dist, '+'))
    }
    d$delRplus[j] <- sample(delRplus.dist, size = 1)
    if(is.na(d$delRplus[j]))
      stop()
    d$Rplus <- cumsum(d$delRplus)
  }
  return(d$Rplus)
}