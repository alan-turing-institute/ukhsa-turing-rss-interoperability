# For reproducibility
set.seed(42)

### Estimate LTLA prevalence ###
library(dplyr)
library(prevdebiasr)
library(parallel)
library(foreach)
library(truncnorm)
source("scripts/SIR_utils.R")

trans_mats <- readRDS("transmats/poisson_SIR_epi_gamma_1.RDS")
vax_df <- readr::read_csv("data/vaccination.csv")
pcr_infectious_df <- readr::read_csv("data/moment_match_infectious.csv")
region_df <- readr::read_csv("data/region.csv") %>%
  left_join(pcr_infectious_df, by = c("phe_region", "mid_week"))
ltla_df <- readr::read_csv("data/ltla.csv") %>%
  left_join(vax_df, by = c("ltla", "mid_week")) %>%
  left_join(select(region_df, -c(Nt, nt, M)),
    by = c("phe_region", "mid_week")
  ) %>%
  select(ltla, phe_region, mid_week, Nt, nt, M, V, alpha, beta)

alpha_testing <- 3e-4
control_debias <- prevdebiasr::get_control_parameters(
    alpha_testing = alpha_testing
)

mid_week_unique <- sort(unique(ltla_df$mid_week))

n_cores <- 2
clust <- makeCluster(n_cores)
doParallel::registerDoParallel(clust)

out_dir <- "output"
control <- c(control_debias, control_SIR)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(control, file.path(out_dir, "control.RDS"), version = 2)

clusterExport(clust, c("control_debias", "control_SIR"))

### PCR positivity, perfect testing ###
imperfect <- FALSE
type <- "PCR_positive"

# Calculate regional bias parameters
delta_df <- region_df %>%
    group_by(phe_region) %>%
    group_modify(~ cbind(
      mid_week = .x$mid_week,
      specify_delta_prior(.x, control_debias, imperfect)
    ))


# Estimate local prevalence for each LTLA
ltla_list <- ltla_df %>%
      left_join(delta_df, by = c("phe_region", "mid_week")) %>%
      group_by(ltla) %>%
      group_split()
ltla_names <- sapply(ltla_list, function(x) x$ltla[1])
ltla_prevalence <- parLapply(clust, ltla_list, local_prevalence,
                            control_debias, imperfect, type)
names(ltla_prevalence) <- ltla_names

# Save output
delta_out_file <- file.path(out_dir, "delta.csv")
readr::write_csv(delta_df, delta_out_file)

ltla_out_file <- file.path(out_dir, "ltla_prevalence.RDS")
saveRDS(ltla_prevalence, ltla_out_file, version = 2)


### Moment-matching for interoperability ###

logit_moments <- list()

for (this_ltla in ltla_names) {
  this_M <- ltla_df %>%
    filter(ltla == this_ltla) %>%
    pull(M) %>%
    first()
  
  I_prop <- control$I_seq / this_M
  
  this_mom <- apply(ltla_prevalence[[this_ltla]]$log_lik, 1, logit_mom, I_prop)
  
  logit_moments[[this_ltla]] <- do.call(rbind, this_mom) %>%
    mutate(ltla = this_ltla, mid_week = mid_week_unique, M = this_M)
}

logit_moment_df <- do.call(bind_rows, logit_moments)

readr::write_csv(logit_moment_df, file.path(out_dir, "/logit_moments.csv"))


# Fit SIR to latest available date
foreach(ltla_name = ltla_names, .packages = c("dplyr", "prevdebiasr")) %dopar% { #ltla_name <- ltla_names[1]#
  source("scripts/SIR_utils.R")
  set.seed(42)
  
  d_ltla <- ltla_df %>%
    filter(ltla == ltla_name)
  I_log_lik <- ltla_prevalence[[ltla_name]]$log_post
  
  d_ltla$V[is.na(d_ltla$V)] <- d_ltla$V[max(which(!is.na(d_ltla$V)))]
  
  SIR_model_out_ltla <- sample_I_R_SIR_from_pre_calc_lik(
    d_ltla, I_log_lik,
    trans_mats, control_debias,
    control_SIR
  )
  
  out_file <- file.path(out_dir, "SIR", paste0(ltla_name, ".RDS"))
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(SIR_model_out_ltla, out_file, version = 2)
}

stopCluster(clust)
