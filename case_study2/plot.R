
library(ggplot2)
library(patchwork)
library(xtable)

# Load outputs of the two models
load('prev_res.RData')
load('prev_res_adj.RData')

# Load LTLA_shp_Reg (contains geography)
load('space_obj.Rdata')

# =========================================================================
# Get parameter estimates (Table 1)
# =========================================================================

#prev_res_adj$summary.fixed
#prev_res$summary.fixed

prev_res_adj$summary.hyperpar
prev_res$summary.hyperpar

xtable::xtable(exp(prev_res_adj$summary.fixed[, c(4,3,5)]) )
xtable::xtable(exp(prev_res$summary.fixed[, c(4,3,5)]) )

# =========================================================================
# Plot (Fig 6 & 7)
# =========================================================================

# Spatial latent field -----------------------------------------------------

latent_field_plot <- function(inla_res, spatial_df, title){
  
  spatial_df$latent_field <- inla_res$summary.random$LTLA_ID$`0.5quant`[1:nrow(LTLA_shp_Reg)]
  
  lat_field_plot <- ggplot(spatial_df) +geom_sf(aes(fill = latent_field), colour = NA) +
    labs(fill = "") + ggtitle(title) + scale_fill_viridis_c() + theme_minimal()
  
  return(lat_field_plot)
}

# Get a plot per model and combine into a single figure  
pp1 <- latent_field_plot(prev_res, LTLA_shp_Reg, 'Naive Model')
pp2 <- latent_field_plot(prev_res_adj, LTLA_shp_Reg, 'Heteroscedastic Model')
combined <- pp1 + pp2 + plot_layout(guides = "collect") & scale_fill_viridis_c()
combined <- combined + plot_annotation(title = 'Latent Spatial Field \U1D709')

ggsave("latent_spatial_field.png")

# Exceedance probabilities (lambda > threshold) --------------------------

exceedance_plot <- function(inla_res, threshold, plot_df, title){
  
  # Retrieve marginals from R-INLA output
  marginals <- inla_res$marginals.random$LTLA_ID
  
  # Get probability marginal is above threshold value
  exc <- sapply(marginals, FUN = function(mar){1-inla.pmarginal(q=threshold, marginal=mar)})
  res <-  cut(exc[1:nrow(LTLA_shp_Reg)],
              breaks = c(0, 0.75, 0.85, 0.95, 1),
              labels = c("[0,0.75)", "[0.75,0.85)", "[0.85,0.95)", "[0.95,1.0]"),
              right=FALSE,
              include.highest=TRUE
  )
  
  # Create plot 
  plot_df$val <-  as.factor(res)
  pp <- ggplot(plot_df) + geom_sf(aes(fill = val), colour = NA) +
    labs(fill = "") + ggtitle(title) +
    scale_fill_brewer(drop=FALSE) + theme_minimal()
  return(pp)
}

# Get a plot per model and combine into a single figure
pp1 <- exceedance_plot(prev_res, 0, LTLA_shp_Reg, 'Naive Model')
pp2 <- exceedance_plot(prev_res_adj, 0, LTLA_shp_Reg, 'Heteroscedastic Model')
combined <- pp1 + pp2 + plot_layout(guides = "collect") & scale_fill_brewer(drop=FALSE)
combined <- combined + plot_annotation(title = 'Exceedance probabilities P(\U1D709 >0)')

ggsave("exceedance.png")


