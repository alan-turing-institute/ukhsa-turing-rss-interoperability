library(dplyr)
source("scripts/plot_utils.R")
source("scripts/SIR_utils.R")

# Run 00_download_data.R and 01_preprocess_data.R first
# ltla_df <- readr::read_csv("data/ltla.csv")

# Link to file in "Interoperability of models" Overleaf, file "shared_data/IR_for_interop.csv"
# IR <- read.csv("https://www.dropbox.com/s/43kjmxq9kqunp22/IR_for_interop.csv?dl=1", stringsAsFactors = F)
# IR <- read.csv("https://www.dropbox.com/s/43kjmxq9kqunp22/Rt-combined.csv?dl=1", stringsAsFactors = F)
IR_w_epimap <- read.csv("https://www.dropbox.com/s/j0gtd2nue8sqfy1/Rt-combined-v2.csv?dl=1", stringsAsFactors = F)

IR_w_epimap$ltla_mid_week <- paste0(IR_w_epimap$ltla, "_", IR_w_epimap$mid_week)
mid_week_unique <- sort(unique(IR_w_epimap$mid_week))
ltla_mid_week_unique <- sort(unique(IR_w_epimap$ltla_mid_week)) 
# IR_w_epimap <- IR_w_epimap[match(ltla_mid_week_unique, IR_w_epimap$ltla_mid_week), ]
epimap_min_date <- as.character(min(IR_w_epimap$Date))#mid_week)
epimap_max_date <- as.character(max(IR_w_epimap$Date))#mid_week)
# str(IR_w_epimap)
# min(IR$mid_week)
# head(IR_w_epimap)
# # Link to file in "Interoperability of models" Overleaf, file "shared_data/Rt_estimates_Epimap_combined.csv"
# rd <- read.csv("https://www.dropbox.com/s/txo5s4wllnzht2o/Rt_estimates_Epimap_combined.csv?dl=1") %>%
#   mutate(date = as.Date(Date)) %>%
#   filter(date <= max(ltla_df$mid_week))
# epimap_min_date <- min(rd$date)
# epimap_max_date <- max(rd$date)

mid_week_unique <- unique(as.character(ltla_df$mid_week))
week_ind <- which(mid_week_unique >= epimap_min_date)
ltla_plot <- sample(ltla_df$ltla, 9)


IR_week <- IR_w_epimap[IR_w_epimap$Date %in% mid_week_unique, ] 
type_curr <- c("epimap_diff_from_others", "debias_diff_from_others", "ED_diff_from_others")[1]
if (type_curr == c("epimap_diff_from_others")) {
  ltla_plot <- unique(IR_week[order(abs(IR_week$epimap_R_m - (IR_week$epimap_debiased_R_m + IR_week$debiased_R_m) / 2), decreasing = T), "ltla"])[5:8]
}
if (type_curr == c("debias_diff_from_others")) {
  ltla_plot <- unique(IR_week[order(abs(IR_week$debiased_R_m - (IR_week$epimap_debiased_R_m + IR_week$epimap_R_m) / 2), decreasing = T), "ltla"])[4:6]
}
if (type_curr == c("ED_diff_from_others")) {
  ltla_plot <- unique(IR_week[order(abs(IR_week$epimap_debiased_R_m - (IR_week$debiased_R_m + IR_week$epimap_R_m) / 2), decreasing = T), "ltla"])[1:3]
}

ltla_unique <- unique(IR$ltla)
set.seed(42)
ltla_plot <- c("Birmingham", "Liverpool", "Havering")

col2rgb("cyan")
dir.create(plot_dir_curr, showWarnings = F)
n_ltla_plot <- length(ltla_plot)
export_plot <- TRUE#TRUE
if (export_plot) {
  pdf(file.path(plot_dir_curr, "interop_R_counts_epimap.pdf"), 12, 7)
}
par(mfcol = c(2, n_ltla_plot + 1), mar = c(1, 1.5, 1, 1), oma = c(3, 0, 2, 4))
pch_ests <- 19
cex_ests <- 1
col_curr <- list(corrected = rgb(0, 0, 1, alpha = .5), uncorrected = rgb(1, 0, 0, alpha = .3))
plot(0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", ty = "n")
trans_alpha <- .25
full_alpha <- .75
full_lwd <- 2
col_use_transparent <- list(debias = rgb(0, 0, 1, alpha = trans_alpha), epimap = rgb(1, 0, 0, alpha = trans_alpha), ED = rgb(0, 1, 1, alpha = trans_alpha))#rgb(.5, 0, .5, alpha = .5))
col_use <- list(debias = rgb(0, 0, 1, alpha = full_alpha), epimap = rgb(1, 0, 0, alpha = full_alpha), ED = rgb(0, 1, 1, alpha = full_alpha))#rgb(.5, 0, .5, alpha = .5))
leg_expand <- 1.3
legend(x = "top", legend = c("De-biasing model", "Epimap model", "Epimap debiased"),
       lty = 1, lwd = 3, col = unlist(col_use), cex = leg_expand, bty = "", title = "Effective R")

plot(0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", ty = "n")
legend(x = "top", legend = c("N = Total", "n = Positive", "n/N = Test positivity"),
       pch = c(19, 1, 3), col = c(1, 1, 4), cex = leg_expand, bty = "", text.col = c(1, 1, 4), title = "Weekly test counts")
cexax = 1
pt_cex <- 1.2
lwd_mid <- 3
ax_expand <- 1.2
for (ltla_curr in ltla_plot) {
  d <- as.data.frame(ltla_df %>%
    filter(ltla == ltla_curr & mid_week >= epimap_min_date & mid_week <= epimap_max_date))
  d$mid_week <- as.character(d$mid_week)
  plot_map <- data.frame(date = seq(as.Date(epimap_min_date), 
                                    as.Date(epimap_max_date), 
                                    by = 1))
  plot_map$day_index <- 1:nrow(plot_map)
  plot_map$date2 <- format(as.Date(plot_map$date), format = "%d/%m/%Y")
  plot_map$date_char <- as.character(plot_map$date)
  xpl <- plot_map$day_index
  
  what_pl <- "R"
  ylim_use <- switch(what_pl, I = c(0, max_prev2), R = c(0, 2.75))
  matplot(x = 0, y = 0, ty = "n", ylim = ylim_use, xlim = range(xpl), 
          xaxt = "n", yaxt = "s", xlab = "", las = 2, ylab = "", xaxs = "i", 
          yaxt = "n", yaxs = "i")
  annotate_months(plot_map, add_axis = F, shade_alpha = .2, for_presentation = T, cex.axis = ax_expand)
  mtext(side = 3, text = substr(ltla_curr, 1, 20), cex = 1.2, line = 1)
  mtext(side = 3, at = xpl[1], text = paste0("(", letters[match(ltla_curr, ltla_plot)], ")"), cex = 1.2, line = 1)
  if(match(ltla_curr, ltla_plot) == 1) {
    mtext(side = 2, line = 3, text = "Effective R", cex = cexax)
    axis(side = 2, cex.axis = ax_expand)
  }
  abline(h = 1)
  vline_at <- NULL
  if(ltla_curr %in% c("Havering", "Redbridge", "Medway", "Thurrock", "Southend-on-Sea", "Waltham Forest")) {
    vline_at <- "2020-12-06"
  }
  if(ltla_curr %in% c("Liverpool")) {
    vline_at <- "2020-11-08"
  }
  abline(v = match(vline_at, plot_map$date_char), lty = 2, lwd = 2)
  xpl_week <- plot_map[match(d$mid_week, plot_map$date_char), "day_index"]
  if (what_pl == "R") {
    matc <- IR_w_epimap[IR_w_epimap$ltla == ltla_curr & IR_w_epimap$mid_week %in% as.character(d$mid_week), c("debiased_R_l", "debiased_R_m", "debiased_R_u")]
    matc <- IR_w_epimap[IR_w_epimap$ltla == ltla_curr & IR_w_epimap$Date %in% plot_map$date_char, c("debiased_R_l", "debiased_R_m", "debiased_R_u")]
  }
  colc <- 1
  lines(plot_map$day_index, matc[, 2], lty = 1, lwd = full_lwd, col = col_use$debias)
  polygon(x = c(plot_map$day_index, rev(plot_map$day_index)), y = c(matc[, 1], rev(matc[, 3])), col = col_use_transparent$debias, border = NA)
  epimap_R <- IR_w_epimap[IR_w_epimap$ltla == ltla_curr & IR_w_epimap$Date %in% plot_map$date_char, c("epimap_R_l", "epimap_R_m", "epimap_R_u", "Date")]
  epimap_R <- epimap_R[match(plot_map$date_char, epimap_R$Date), ]
  lines(plot_map$day_index, epimap_R[, 2], lty = 1, lwd = full_lwd, col = col_use$epimap)
  polygon(x = c(plot_map$day_index, rev(plot_map$day_index)), y = c(epimap_R[, 1], rev(epimap_R[, 3])), col = col_use_transparent$epimap, border = NA)
  ED_R <- IR_w_epimap[IR_w_epimap$ltla == ltla_curr & IR_w_epimap$Date %in% plot_map$date_char, c("epimap_debiased_R_l", "epimap_debiased_R_m", "epimap_debiased_R_u", "Date")]
  ED_R <- ED_R[match(plot_map$date_char, ED_R$Date), ]
  lines(plot_map$day_index, ED_R[, 2], lty = 1, lwd = full_lwd, col = col_use$ED)
  polygon(x = c(plot_map$day_index, rev(plot_map$day_index)), y = c(ED_R[, 1], rev(ED_R[, 3])), col = col_use_transparent$ED, border = NA)
  mid_week_plot <- mid_week_unique[mid_week_unique %in% IR_w_epimap$Date]
  ypl <- d[match(mid_week_plot, d$mid_week), c("Nt", "nt")]
  matplot(match(mid_week_plot, plot_map$date_char), ypl, ty = "b",
          log = "", col = 1, lwd = c(1, 1), lty = 1, xaxt = "n", pch = c(19, 1), ylab = "", las = 2, 
          yaxt = "n", cex = pt_cex,
          ylim = c(0, max(ypl)))#,
  axis(side = 2, cex.axis = ax_expand)
  if(match(ltla_curr, ltla_plot) == 1) {
    mtext(side = 2, line = 3, text = "Test counts", cex = cexax)
  }
  annotate_months(plot_map, add_axis = T, shade_alpha = .2, for_presentation = T, cex.axis = ax_expand)
  par(new = T)
  d$case_pos <- d$nt / d$Nt
  plot(match(mid_week_plot, plot_map$date_char), d[match(mid_week_plot, d$mid_week), c("case_pos")], ty = "b",
          log = "", lty = 1, xaxt = "n", yaxt = "n", col = "blue", ylim = c(0, .5), pch = 3, ylab = "")
  abline(v = match(vline_at, plot_map$date_char), lty = 2, lwd = 2)
  if(match(ltla_curr, ltla_plot) == n_ltla_plot) {
    mtext(side = 4, line = 3,  text = "Test positivity", col = 4, cex = cexax)
    axis(side = 4, col.axis = "blue", col = "blue", cex.axis = ax_expand)
  }
}
if (export_plot) {
  dev.off()
}

map_look <- IR_w_epimap[IR_w_epimap$Date == "2020-12-07", ]
map_look[order(map_look$epimap_debiased_R_m), ]
north_warks_R_E <- map_look[map_look$ltla == "North Warwickshire", "epimap_R_m"]
north_warks_R_ED <- map_look[map_look$ltla == "North Warwickshire", "epimap_debiased_R_m"]
north_warks_R_D <- map_look[map_look$ltla == "North Warwickshire", "debiased_R_m"]
craven_R_E <- map_look[map_look$ltla == "Craven", "epimap_R_m"]
craven_R_ED <- map_look[map_look$ltla == "Craven", "epimap_debiased_R_m"]
craven_R_D <- map_look[map_look$ltla == "Craven", "debiased_R_m"]

liverpool <- ltla_df[ltla_df$ltla == "Liverpool", ]
liverpool$test_pos <- liverpool$nt / liverpool$Nt 
# dir_text_numbers <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/text_numbers/cut_vs_full_comp"
tes_pos_liv_week_1 <- liverpool$test_pos[liverpool$mid_week == "2020-11-01"] * 100
tes_pos_liv_week_2 <- liverpool$test_pos[liverpool$mid_week == "2020-11-08"] * 100

# dir_text_numbers <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Interoperability of models/text_numbers/epimap_case_study"
# dir.create(dir_text_numbers, showWarnings = F)

save_num <- c("tes_pos_liv_week_1", "tes_pos_liv_week_2")
for(numc in save_num)
  write.table(formatC(eval(as.name(numc)), format = "f", digits = 1), file = paste(dir_text_numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)


save_num <- c("north_warks_R_E", "north_warks_R_ED", "north_warks_R_D", "craven_R_E", "craven_R_ED", "craven_R_D")
for(numc in save_num)
  write.table(formatC(eval(as.name(numc)), format = "f", digits = 2), file = paste(dir_text_numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)






