

IR$ltla_week <- paste0(IR$ltla, "_", IR$mid_week)
IR[, paste0("Ix_", c("l", "m", "u"))] <- NA
# Ix_all <- data.frame()
for (ltla_curr in ltla_unique) {
  M_curr <- ltla_df[match(ltla_curr, ltla_df$ltla), "M"]
  add <- data.frame(ltla = ltla_curr, mid_week = mid_week_unique)
  indices_quants_curr_raw <- t(apply(ltla_prevalence[[ltla_curr]]$norm_post, 1, function(v) findInterval(c(.025, .5, .975), cumsum(v))))
  indices_quants_curr <- indices_quants_curr_raw
  indices_quants_curr[, 1] <- pmax(1, indices_quants_curr_raw[, 1] - 1)
  indices_quants_curr[, 3] <- pmin(length(control$I_seq), indices_quants_curr_raw[, 3] + 1)
  # pi_quants_curr <- indices_quants_curr
  # pi_quants_curr[] <- control$I_seq[indices_quants_curr] / M_curr * 100
  IR[match(paste0(ltla_curr, "_", mid_week_unique), IR$ltla_week), paste0("Ix_", c("l", "m", "u"))] <- 
    control$I_seq[indices_quants_curr] / unlist(M_curr) * 100#pi_quants_curr
  # names(I_quants_curr) <- paste0(c("l", "m", "u"), "_Ix")
}

str(ltla_prevalence[[ltla_curr]])

phe_region_unique
# ltlas_focus <- unlist(unique(ltla_df[ltla_df$phe_region == "London", "ltla"]))
ltlas_focus <- unlist(unique(ltla_df[match(phe_region_unique[c(1, 2, 5)], ltla_df$phe_region), "ltla"]))
ltla_curr <- ltlas_focus[1]
irc <- IR[IR$ltla == ltla_curr, ]
panel_count <- 1
max_prev_I <- max(IR$I_u[IR$ltla %in% ltlas_focus])

plot_map <- data.frame(date = seq(min(as.Date(ltla_df$mid_week)) - 3, 
                                  max(as.Date(ltla_df$mid_week)) + 3, 
                                  by = 1))
plot_map$day_index <- 1:nrow(plot_map)
xpl <- plot_map$day_index
ltla_curve <- ltlas_focus#[1:4]
col_reg <- rainbow(length(ltla_curve), alpha = .75)
names(col_reg) <- ltla_curve
# plot_dir <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/figures"
jpeg(paste0(plot_dir, "/SIR_interoperability.jpeg"), 9, 6, res = 750, units = "in")
par(mfcol = c(2, length(ltla_curve)), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(6, 5, 6, 5))
for (ltla_curr in ltla_curve) {
  for (what_pl in c("I", "R")[1:2]) {
    # panel_count <- panel_count + 1
    ylim_use <- switch(what_pl, I = c(0, max_prev_I), R = c(0, 2.5))
    matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
            xaxt = "n", yaxt = "s", xlab = "", las = 2, ylab = "", xaxs = "i", yaxt = "n")
    if(ltla_curr == ltla_curve[1])
      axis(2, las = 1, cex.axis = 1.1)
      # axis(ifelse(what_pl == "I", 2, 4), las = 1)
    if(what_pl == "R")
      abline(h = 1)
  
    if(ltla_curr == ltla_curve[1]) {
      mtext(side = 2, text = switch(what_pl,
                                  I = "% Prevalence",
                                  R = "Effective R"), line = 3)
    }
    if(what_pl == "I") {
      matc <- IR[IR$ltla == ltla_curr, c("I_l", "I_m", "I_u")]
      matc_x <- IR[IR$ltla == ltla_curr, c("Ix_l", "Ix_m", "Ix_u")]
    }
    if(what_pl == "R") {
      matc <- IR[IR$ltla == ltla_curr, c("R_l", "R_m", "R_u")]
    }
    
    xpl_week <- plot_map[match(mid_week_unique, plot_map$date), "day_index"]
    colc <- rgb(0, 0, 1, alpha = .5)#"grey"#col_reg[ltla_curr]
    if(what_pl == "I") {
      points(xpl_week, matc_x[, 2], pch = 19, col = 1)
      for (i in 1:length(xpl_week)) {
        lines(rep(xpl_week[i], 2), matc_x[i, c(1, 3)], lty = 1, col = 1)
      }
    }
    lines(xpl_week, matc[, 2], lty = 1, lwd = 2, col = colc)
    lines(xpl_week, matc[, 3], lty = 1, lwd = 1, col = colc)
    lines(xpl_week, matc[, 1], lty = 1, lwd = 1, col = colc)
    annotate_months(plot_map, add_axis = (what_pl == "R"), shade_alpha = .2, for_presentation = T)
    if(what_pl == "I")
      mtext(side = 3, text = ltla_curr, cex = 1, line = 1)
  }
}
dev.off()

