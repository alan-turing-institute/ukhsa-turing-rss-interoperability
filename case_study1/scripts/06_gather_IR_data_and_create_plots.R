library(dplyr)
library(rgdal)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

options(bitmapType = "cairo-png")

source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

id <- "PCR_positive_Perfect"
control_par_id <- "AR0.99sd1Rsd0.2"
output_dir <- "output"
out_dir <- file.path(output_dir, control_par_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_data_from_archive <- TRUE
data_dir_to_use <- ifelse(load_data_from_archive, "data_from_archive", "data")
out_dir_to_use <- ifelse(load_data_from_archive, file.path(data_dir_to_use, "output", control_par_id),
                         file.path("output", control_par_id))

output_to_overleaf <- F#TRUE
switching_index <- ifelse(output_to_overleaf, 2, 1)
plot_dir_curr <- c("plots", "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/figures")[switching_index]
dir_text_numbers <- c(file.path("output", "numbers_for_text"),
                      "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/text_numbers/cut_vs_full_comp")[switching_index]
dir.create(plot_dir_curr, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_text_numbers, recursive = TRUE, showWarnings = FALSE)

##########################
### Load data & output ###
##########################
ltla_df <- readr::read_csv(file.path(data_dir_to_use, "ltla.csv"))
react_date_df <- readr::read_csv(file.path(data_dir_to_use, "react_dates.csv"))
react_ltla_df <- readr::read_csv(file.path(data_dir_to_use, "react_ltla.csv"))
raw_pillar2_df <- readr::read_csv(file.path(data_dir_to_use, "raw_pillar2.csv"))

delta_df <- readr::read_csv(file.path(out_dir_to_use, id, "delta.csv"))
ltla_prevalence <- readRDS(file.path(out_dir_to_use, "ltla_prevalence_pcr_perfect.RDS"))
ltla_pop <- ltla_df %>%
  distinct(ltla, M)

# LTLA_shp_Reg <- get_ltla_shape_file()

control <- prevdebiasr::get_control_parameters()
# Quantiles to plot
quant_plot <- c(0.025, 0.5, 0.975)
react_rounds <- sort(unique(react_ltla_df$round))


organise_new_SIR_output <- FALSE


ltla_unique <- unique(ltla_df$ltla)
mid_week_unique <- unique(ltla_df$mid_week)
n_weeks <- length(mid_week_unique)
all_mid_week_cut <- mid_week_unique
date_recent <- max(mid_week_unique)

if (!organise_new_SIR_output) {
  IR <- readr::read_csv(file.path("data_from_archive", "IR_for_interop.csv"))
  ltla_unique <- unique(IR$ltla)
}
if (organise_new_SIR_output) {
  # Load SIR output
  out_files <- list.files(file.path(out_dir_to_use, id, "SIR"), full.names = TRUE)
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
}


#################################
### Debiased prevalence plots ###
#################################
source("scripts/06a_plot_cut_vs_full.R")
source("scripts/06b_plot_SIR.R")
source("scripts/06c_plot_epimap_comparison.R")

