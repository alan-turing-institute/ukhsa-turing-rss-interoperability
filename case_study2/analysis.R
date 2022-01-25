library(dplyr)
library(tidyr)
library(INLA)

# =========================================================================
# Load and format data
# =========================================================================

#source("data_prep.R")

## File produced by running `data_prep.R`
prev_fin <- readRDS(file='prev_data.rds')

prev_fin <- prev_fin[!is.na(prev_fin$week),]

prev_fin$date_LTLA_ID <- as.numeric(interaction(prev_fin$LTLA_ID, prev_fin$date_ID))

# Get standard deviations by LTLA for heteroscedastic model
prev_sd <- prev_fin %>% dplyr::select(date_LTLA_ID, sd) %>% distinct()
prev_sd <- prev_sd[order(prev_sd$date_LTLA_ID), ]

# =========================================================================
# Models
# =========================================================================

# Naive model ------------------------------------------------------------

formula <- mean ~ 1  + IMD_stand + bame_stand + 
  f(date_ID, model = "rw2") +
  f(LTLA_ID, model = "bym2", graph = "data/W.adj", scale.model = TRUE, constr =
      TRUE, hyper = list(theta1 = list("PCprior", c(1,   0.01)), 
                         theta2 = list("PCprior", c(0.5, 0.5))))


prev_res <- inla(formula, data = prev_fin, family = "gaussian",
                  verbose = TRUE, 
                  num.threads = 2,
                  control.fixed=list(prec=1,prec.intercept=1),
                  control.inla=list(int.strategy="eb", strategy="adaptive"))

summary(prev_res)
save(prev_res, file="prev_res.RData")

# Heteroscedastic model ----------------------------------------------------

hyper_prior <- list(prec = list(prior = "gamma", initial = 1, fixed = TRUE))

formula_adj <- mean ~ 1  + IMD_stand + bame_stand + 
  f(date_ID, model = "rw2") +
  f(LTLA_ID, model = "bym2", graph = "data/W.adj", scale.model = TRUE, constr =
      TRUE, hyper = list(theta1 = list("PCprior", c(1,   0.01)), theta2 =
                           list("PCprior", c(0.5, 0.5)))) +
  f(date_LTLA_ID, model = "generic0", 
    Cmatrix = Diagonal(nrow(prev_sd), 1/prev_sd$sd^2), hyper = hyper_prior ) 


prev_res_adj <- inla(formula_adj, data = prev_fin, family = "gaussian",
                     verbose = TRUE, 
                     num.threads = 2,
                     control.fixed=list(prec=1,prec.intercept=1),
                     control.inla=list(int.strategy="eb", strategy="adaptive"))
summary(prev_res_adj)
save(prev_res_adj, file="prev_res_adj.RData")

