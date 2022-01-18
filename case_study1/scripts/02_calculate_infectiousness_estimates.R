library(dplyr)

region_df <- readr::read_csv("data/region.csv")

#=========================================================================
# Load posterior samples for P(PCR+ | time since infection)
# Samples are from Fig 3b in https://www.medrxiv.org/content/10.1101/2020.11.24.20229948v3.full.pdf
# Code to produce posterior samples: https://github.com/cmmid/pcr-profile
# 
# NOTES:
# The model is fit to days 0-30 (inclusive) in increments of 0.1
# This means we have 301 columns (the time index) and 4000 posterior samples (rows) for each
#=========================================================================

pcr_df <- readr::read_csv("data/samples_pcr.csv")
pcr_samples <- pcr_df[,2:302]

#=========================================================================
# p(infectious | infected) && p(pcr+ | infected)
#=========================================================================
# NOTE: here we index week k 1-4 rather than 0-3 as in the paper

# 1: p(infectious | infected) for each week k, eq.40
multipliers = c(6/7, 5/7, 0, 0)

# 2: p(pcr+ | infected) for each week k, eq.42
get_pcr_infect_week <- function(k, pcr_samples){
  
  start_week <- 1 + (k-1) * 70
  end_week <- start_week + 70
  mean_samples <- rowMeans(pcr_samples[, start_week:end_week])
  
  return(mean_samples)
}

#=========================================================================
# P(Infectious | PCR+), constant incidence
#=========================================================================
#dists <- pcr_samples / rowSums(pcr_samples)
#pcrplus_dist <- rowSums(dists[,11:110])

for (k in 1:4) {
  
  mean_samples <- get_pcr_infect_week(k, pcr_samples)
  
  # eq. 39
  denom = mean_samples 
  num = denom * multipliers[k]
  
  if (k==1){
    nums <- num
    denoms <- denom
  } else{
    nums <- nums + num
    denoms <- denoms + denom
  }
  
}

pcrplus_dist <- nums/denoms

df <- data.frame(samples = pcrplus_dist)
mu <- mean(df$samples)
sd <- sd(df$samples)
low_ci <- quantile(df$samples, prob=0.025)
up_ci <- quantile(df$samples, prob=0.975)

#=========================================================================
# P(Infectious | PCR+), adjusted for incidence
# NEW METHOD
#=========================================================================

# p(infected) depends on incidence and is calculated for each region & week combination

dates <- unique(region_df$mid_week)
phe_regions <- unique(region_df$phe_region)

# region loop
df_list <- list()
for (r in phe_regions) {
  
  r_df <- region_df %>% 
    filter(phe_region == r) %>% 
    arrange(mid_week)
  counts <- r_df$nt
  n_counts = length(counts)
  
  nums <- 0
  denoms <- 0
  
  # lag loop
  for (k in 1:4){
    
    mean_samples <- get_pcr_infect_week(k, pcr_samples)
    
    # 3: p(infected), eq. 45
    inc_ests = c()
    # loop through mid_week dates
    for (i in 4:n_counts){
      start = i - 3
      curr_slice = counts[start:i]
      roll = sum(counts[start:i])
      idx = length(curr_slice) - k + 1
      inc_ests <- c(inc_ests, curr_slice[idx] / roll)
    }
    
    # eq. 39
    denom = outer(mean_samples, inc_ests)
    num = denom * multipliers[k]
    
    if (k==1){
      nums <- num
      denoms <- denom
    } else{
      nums <- nums + num
      denoms <- denoms + denom
    }
    
  }
  result = nums / denoms
  
  res_df <- data.frame(
    upper = c(rep(up_ci, 3), apply(result, 2, quantile, prob = 0.975, na.rm = T)),
    lower = c(rep(low_ci, 3), apply(result, 2, quantile, prob = 0.025, na.rm = T)),
    mean = c(rep(mu, 3), apply(result, 2, mean, na.rm = T)),
    var = c(rep(sd ^ 2, 3), apply(result, 2, var, na.rm = T)),
    mid_week = dates
  )
  
  res_df$phe_region <- r
  res_df$alpha <- (res_df$mean ^ 2 * (1 - res_df$mean)) / res_df$var - res_df$mean
  res_df$beta <- (res_df$mean * (1 - res_df$mean) / res_df$var - 1) * (1 - res_df$mean)
  df_list[[r]] <- res_df
}

results_df <- do.call(rbind, df_list)
rownames(results_df) <- NULL
readr::write_csv(results_df, "data/moment_match_infectious.csv")

