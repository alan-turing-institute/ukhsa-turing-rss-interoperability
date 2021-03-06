library(dplyr)
library(tidyr)

source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")


out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (Sys.info()["user"] == "nicho") {
  plot_dir <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/figures"
  dir_text_numbers_epimap <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/text_numbers/epimap_case_study"
  dir_text_numbers_case_study1 <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/text_numbers/cut_vs_full_comp"
} else {
  plot_dir <- "plots"
  dir_text_numbers_epimap <- paste0(out_dir, "/text_numbers/epimap_case_study")
  dir_text_numbers_case_study1 <- paste0(out_dir, "/text_numbers/cut_vs_full_comp")
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_text_numbers_epimap, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_text_numbers_case_study1, showWarnings = FALSE, recursive = TRUE)


##########################
### Load data & output ###
##########################
ltla_df <- readr::read_csv("data/ltla.csv")
react_date_df <- readr::read_csv("data/react_dates.csv")
react_ltla_df <- readr::read_csv("data/react_ltla.csv")
raw_pillar2_df <- readr::read_csv("data/raw_pillar2.csv")

delta_df <- readr::read_csv(file.path(out_dir, "delta.csv"))
ltla_prevalence <- readRDS(file.path(out_dir, "ltla_prevalence.RDS"))
ltla_pop <- ltla_df %>%
  distinct(ltla, M)

control <- prevdebiasr::get_control_parameters()

quant_plot <- c(0.025, 0.5, 0.975)
react_rounds <- sort(unique(react_ltla_df$round))

###########################
# Load SIR output
ltla_unique <- unique(ltla_df$ltla)
mid_week_unique <- unique(ltla_df$mid_week)
n_weeks <- length(mid_week_unique)
all_mid_week_cut <- mid_week_unique
date_recent <- max(mid_week_unique)

out_files <- list.files(file.path(out_dir, "SIR"), full.names = TRUE)
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
readr::write_csv(IR, file.path(out_dir, "IR_for_interop.csv"))


source("scripts/06a_cut_vs_full.R")
source("scripts/06b_SIR_plot.R")
source("scripts/06c_epimap_comparison_plots.R")
