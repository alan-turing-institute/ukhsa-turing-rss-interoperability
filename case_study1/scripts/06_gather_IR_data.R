library(dplyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

options(bitmapType = "cairo-png")

source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

id <- "PCR_positive_perfect"
control_par_id <- "AR0.99sd1Rsd0.2"
# out_dir <- file.path("output", control_par_id)
# plot_dir <- file.path(out_dir, id, "plots")
plot_dir <- "plots"
output_dir <- "output"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

##########################
### Load data & output ###
##########################
ltla_df <- readr::read_csv("data/ltla.csv")
react_date_df <- readr::read_csv("data/react_dates.csv")
react_ltla_df <- readr::read_csv("data/react_ltla.csv")
raw_pillar2_df <- readr::read_csv("data/raw_pillar2.csv")

delta_df <- readr::read_csv(file.path(out_dir, "delta_pcr_perfect.csv"))
ltla_prevalence <- readRDS(file.path(out_dir, "ltla_prevalence_pcr_perfect.RDS"))
ltla_pop <- ltla_df %>%
  distinct(ltla, M)

LTLA_shp_Reg <- get_ltla_shape_file()

control <- prevdebiasr::get_control_parameters()

# Quantiles to plot
quant_plot <- c(0.025, 0.5, 0.975)
react_rounds <- sort(unique(react_ltla_df$round))

# Load SIR output
ltla_unique <- unique(ltla_df$ltla)
mid_week_unique <- unique(ltla_df$mid_week)
n_weeks <- length(mid_week_unique)
all_mid_week_cut <- mid_week_unique
date_recent <- max(mid_week_unique)

out_files <- list.files(file.path(out_dir, id, "SIR"), full.names = TRUE)
SIR_model_results <- lapply(out_files, readRDS)
names(SIR_model_results) <- sub(".RDS", "", basename(out_files))

IR <- data.frame()
Rl <- Il <- list()

for (ltla_curr in ltla_unique) {
  
  this_M <- ltla_pop %>%
    filter(ltla == ltla_curr) %>%
    pull(M)
  
  saml_biased <- SIR_model_results[[ltla_curr]]
  its_keep <- (control_SIR$burn_in + 1):control_SIR$n_iters
  I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, 
                          function(v) quantile(v, qplot, na.rm = T))) / this_M * 100
  SIR_model_results[[ltla_curr]]$I_quant <- I_quant_curr
  SIR_model_results[[ltla_curr]]$R_quant <- t(apply(saml_biased$R[, its_keep], 1, 
                                                    function(v) quantile(v, qplot, na.rm = T)))
  SIR_model_results[[ltla_curr]]$R_quant[n_weeks, ]
  R_add <- as.data.frame(SIR_model_results[[ltla_curr]]$R_quant)
  I_add <- as.data.frame(SIR_model_results[[ltla_curr]]$I_quant)
  rownames(R_add)
  names(R_add) <- paste0("R_", c("l", "m", "u"))
  names(I_add) <- paste0("I_", c("l", "m", "u"))
  add_all <- cbind(I_add, R_add)
  add_all$namdat <- rownames(add_all)
  add_all$ltla <- ltla_curr
  #add_all$mid_week <- sapply(strsplit(add_all$namdat, "_"), function(x) x[2])
  add_all$mid_week <- mid_week_unique
  IR <- rbind(IR, add_all)
}
R_all <- IR %>%
  filter(mid_week == date_recent) %>%
  select(ltla, l = R_l, m = R_m, u = R_u)
I_all <- IR %>%
  filter(mid_week == date_recent) %>%
  select(ltla, l = I_l, m = I_m, u = I_u)
rownames(R_all) <- R_all$ltla
rownames(I_all) <- I_all$ltla

str(IR)
readr::write_csv(IR, path = file.path(output_dir, "IR_for_interop.csv"))

