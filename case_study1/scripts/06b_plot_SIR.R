#####################################################
# Run 06_gather_IR_data first!!
IR <- as.data.frame(IR)
IR$ltla_week <- paste0(IR$ltla, "_", IR$mid_week)
IR[, paste0("Ix_", c("l", "m", "u"))] <- NA
for (ltla_curr in ltla_unique) {
  M_curr <- ltla_df[match(ltla_curr, ltla_df$ltla), "M"]
  add <- data.frame(ltla = ltla_curr, mid_week = mid_week_unique)
  indices_quants_curr_raw <- t(apply(ltla_prevalence[[ltla_curr]]$norm_post, 1, function(v) findInterval(c(.025, .5, .975), cumsum(v))))
  indices_quants_curr <- indices_quants_curr_raw
  indices_quants_curr[, 1] <- pmax(1, indices_quants_curr_raw[, 1] - 1)
  indices_quants_curr[, 3] <- pmin(length(control$I_seq), indices_quants_curr_raw[, 3] + 1)
  IR[match(paste0(ltla_curr, "_", mid_week_unique), IR$ltla_week), paste0("Ix_", c("l", "m", "u"))] <- 
    control$I_seq[indices_quants_curr] / unlist(M_curr) * 100
}

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
ltla_curve <- ltlas_focus
col_reg <- rainbow(length(ltla_curve), alpha = .75)
names(col_reg) <- ltla_curve
pdf(paste0(plot_dir_curr, "/SIR_interoperability.pdf"), 9, 6)#, res = 750, units = "in")
# pdf("C:/Temp/test.pdf")
par(mfcol = c(2, length(ltla_curve)), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(6, 5, 6, 5))
for (ltla_curr in ltla_curve) {
  for (what_pl in c("I", "R")[1:2]) {
    ylim_use <- switch(what_pl, I = c(0, max_prev_I), R = c(0, 2.5))
    matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
            xaxt = "n", yaxt = "s", xlab = "", las = 2, ylab = "", xaxs = "i", yaxt = "n")
    if(ltla_curr == ltla_curve[1])
      axis(2, las = 1, cex.axis = 1.1)
    if(what_pl == "R")
      abline(h = 1)
    if(ltla_curr == ltla_curve[1]) {
      mtext(side = 2, text = switch(what_pl,
                                  I = "% Prevalence",
                                  R = "Effective R"), line = 3)
    }
    d <- ltla_df %>%
      filter(ltla == ltla_curr)
    xpl_week <- plot_map[match(d$mid_week, plot_map$date), "day_index"]
    if(what_pl == "I") {
      matc <- IR[IR$ltla == ltla_curr, c("I_l", "I_m", "I_u")]
      matc_x <- IR[IR$ltla == ltla_curr, c("Ix_l", "Ix_m", "Ix_u")]
    }
    if(what_pl == "R") {
      matc <- IR[IR$ltla == ltla_curr, c("R_l", "R_m", "R_u")]
    }
    colc <- rgb(0, 0, 1, alpha = .5)
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

# plot_map
# pdf("C:/Temp/test.pdf")
# matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
#         xaxt = "n", yaxt = "s", xlab = "", las = 2, ylab = "", xaxs = "i", yaxt = "n")
# 
# d <- plot_map
# shade_alpha <- .5
# d$month <- months(d$date)
# d$month <- format(d$date, "%Y-%m")
# month_unique <- unique(d$month)  
# month_name <- format(strptime(paste0(month_unique, "-01"), format = "%Y-%m-%d"), "%b")
# names(month_name) <- month_unique
# greys_shade <- grey(c(0.75, 1), alpha = shade_alpha)
# month_shade_col <- greys_shade[1:length(month_unique) %% 2 + 1]
# names(month_shade_col) <- month_unique
# month_mean_index <- sapply(month_unique, function(m) mean(which(d$month == m)))
# month_lower_divider <- sapply(month_unique, function(m) min(which(d$month == m)) - .5)
# month_upper_divider <- sapply(month_unique, function(m) max(which(d$month == m)) + .5)
# for(month_shade in month_unique){
#   x_coords <- c(rep(month_lower_divider[month_shade], 2), rep(month_upper_divider[month_shade], 2), month_lower_divider[month_shade])
#   big_num <- 10
#   y_coords <- big_num * c(-1, 1, 1, -1, -1)
#   polygon(x = x_coords, y = y_coords, col = month_shade_col[month_shade], border = NA)  
# }
# dev.off()
# if (add_axis) {
#   labs <- month_name[month_unique]
#   ats <- month_mean_index
#   #   if(for_presentation) {
#   #      labs[labs == "May"] <- ""
#   #    }
#   axis(side = 1, at = ats, labels = labs, las = 2, tcl = 0, ...)
# }