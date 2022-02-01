library(patchwork)
library(dplyr)
library(readxl)
library(readr)
library(lubridate)
library(RColorBrewer)
options(bitmapType = "cairo-png", device = "X11")

## Local Authority Districts (same as LTLA, checked) to Regions lookup:
LTLA2Reg <- read_csv("data/Local_Authority_District_to_Region__December_2020__Lookup_in_England.csv")

out_files <- list.files("output/SIR", full.names=TRUE)
SIR_model_results <- lapply(out_files, readRDS)
ltla_df <- readr::read_csv("data/ltla.csv")
ltla_unique <- unique(ltla_df$ltla)
mid_week_unique <- unique(ltla_df$mid_week)
n_weeks <- length(mid_week_unique)
all_mid_week_cut <- mid_week_unique
date_recent <- max(mid_week_unique)
path_to_ons_pop <- "data/ukmidyearestimates20192020ladcodes.xls"
pop_size <- read_excel(path_to_ons_pop, sheet = "MYE2 - Persons", skip = 3) %>%
  dplyr::select(1:4) %>%
  filter(Name %in% ltla_unique)
#  filter(Code %in% c(voc_df$LAD20CD, "E06000060"))
northamptonshire_pop <- readr::read_csv("data/northamptonshire_pop_nomis.csv", skip = 5) %>%
  rename(ltla = 1, M = 2)
ltla_pop <- pop_size %>%
  dplyr::select(ltla = Name, M = `All ages`) %>%
  filter(ltla != "Northamptonshire") %>%
  bind_rows(northamptonshire_pop)

voc_df <- readr::read_csv("data/variants.csv") %>%
  rename(denom = n_tot) %>%
  filter(denom >= 10) %>%
  left_join(select(LTLA2Reg, LAD20NM, LAD20CD), by = c("ltla" = "LAD20NM")) %>%
  mutate(week = week(mid_week),
         `Percent delta VOC` = 100 * prop_delta) %>%
  filter(!is.na(ltla))
### Load prevalence output
names(SIR_model_results) <- sub(".RDS", "", basename(out_files))
IR <- data.frame()
R_all <- I_all <- data.frame()
Rl <- Il <- list()
control <- list(n_iters = 200, burn_in = 100, qplot = c(0.025, 0.5, 0.975))
for (ltla_curr in ltla_unique) {
  this_M <- ltla_pop %>%
    filter(ltla == ltla_curr) %>%
    pull(M)
  saml_biased <- SIR_model_results[[ltla_curr]]
  its_keep <- (control$burn_in + 1):control$n_iters
  I_quant_curr <- t(apply(saml_biased$I[, its_keep], 1, 
                          function(v) quantile(v, control$qplot, na.rm = T))) / this_M * 100
  SIR_model_results[[ltla_curr]]$I_quant <- I_quant_curr
  SIR_model_results[[ltla_curr]]$R_quant <- t(apply(saml_biased$R[, its_keep], 1, 
                                                    function(v) quantile(v, control$qplot, na.rm = T)))
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

### PLOTS ###

# Plotting libraries 
library(rgdal); library(sf); library(dplyr); library(ggplot2); library(gridExtra); library(grid); library(lattice)
#install.packages("ggdraw")
# import LTLA shapefile
LTLA_shp = readOGR("data/Local_Authority_Districts_(May_2020)_Boundaries_UK_BFE-shp/Local_Authority_Districts_(May_2020)_Boundaries_UK_BFE.shp")
LTLA_shp <- st_as_sf(LTLA_shp)
# remove LTLAs for which we don't have testing data (i.e. all that are not in England)
LTLA_shp_EN = LTLA_shp %>% 
  dplyr::filter(lad20cd %in% LTLA2Reg$LAD20CD)
# add region lookup to LTLA shapefile
LTLA_shp_Reg  <- left_join(LTLA_shp_EN, LTLA2Reg, by = c("lad20cd" = "LAD20CD"))
bar_width <- 0.35 # cm
xmin <- 82672; xmax <- 655653.8
ymin <- 5342.7; ymax <- 657536
yoffset <- (1/5)*(ymax-xmin); xoffset <- (1/5)*(xmax-xmin)

mid_week_start <- as.Date("2021-04-18")
mid_weeks <- mid_week_start + (1:5 * 7)
mid_weeks <- mid_week_start + (0:3 * 21)
mid_weeks <- mid_week_start - 7 + (1:4 * 14)

var_plots <- vector("list", 5)
for (i in seq_along(mid_weeks)) {
  
  this_mid_week <- mid_weeks[i]
  variant_df_week <- dplyr::filter(voc_df, week == week(this_mid_week))
  var_week_shp <- left_join(LTLA_shp_Reg, variant_df_week, by = c("lad20cd" = "LAD20CD"))
  var_plt <- ggplot(data = var_week_shp) + 
    geom_sf(aes(fill = `Percent delta VOC`), color = NA) + 
    scale_fill_viridis_c(limits=c(0,100), 
                         guide = guide_colourbar(title.position = "top"),
                         na.value = "gray90") +
    theme_void()
  
  var_lnd <- ggplot(data = filter(var_week_shp, RGN20NM == "London")) + 
    geom_sf(aes(fill = `Percent delta VOC`), color = NA, show.legend = FALSE) + 
    scale_fill_viridis_c(limits = c(0, 100),
                         na.value = "gray90") +
    theme_void() + 
    ggtitle("London") + 
    theme(plot.title = element_text(size = 5))
  
  var_plots[[i]] <- var_plt +
    annotation_custom(grob = ggplotGrob(var_lnd), xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                      ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) + 
    theme_void() +
    theme(plot.title = element_text(size = 5), 
          legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.key.height = unit(0.2, 'cm'), #change legend key height
          legend.key.width = unit(bar_width, 'cm'), #change legend key width
          legend.title = element_text(size=7), #change legend title font size
          legend.text = element_text(size=5)) #
}

n_week_plot <- 20
maxy_plot <- 2
rng <- c(0, maxy_plot)
#rev(rev(mid_week_unique)[1 + 1 * (0:(n_week_plot - 1))])
high_col = "red"
mid_col = "yellow"
low_col = "green"
midpoint = 1
breaks = seq(rng[1], rng[2], .5)
leg_name = "Rt"
ltla_R_plots <- vector("list", 5)
for (i in seq_along(mid_weeks)) {
  week_plot_curr = mid_weeks[i]
  IR_week <- IR[IR$mid_week == week_plot_curr, ]
  curr_y <- IR_week[match(LTLA_shp_Reg$LAD20NM, IR_week$ltla), "R_m"]
  curr_y[curr_y > maxy_plot] <- maxy_plot
  curr_y[is.na(curr_y)] <- mean(curr_y, na.rm = T)
  #LTLA_shp_Reg2 <- LTLA_shp_Reg
  #LTLA_shp_Reg2$lf_mean <- curr_y
  #LTLA_shp_Reg2$lf_mean[56] <- mean(LTLA_shp_Reg2[fill_gaps[[56]],]$lf_mean, na.rm = TRUE)
  LTLA_shp_Reg$lf_mean <- curr_y
  y_lnd <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "London")) + 
    geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE) + 
    theme_void() + ggtitle("London") + 
    scale_fill_distiller(palette = "RdYlBu", limits =c(0, 2)) + theme(plot.title = element_text(size = 5), 
                                                                      legend.key.size = unit(0.2, 'cm'), #change legend key size
                                                                      legend.key.height = unit(0.2, 'cm'), #change legend key height
                                                                      legend.key.width = unit(0.2, 'cm'), #change legend key width
                                                                      legend.title = element_text(size=7), #change legend title font size
                                                                      legend.text = element_text(size=5)) #change legend text font size
  lnd_grob <- ggplotGrob(y_lnd)
  ltla_R_plots[[i]] <-  ggplot(data = LTLA_shp_Reg) + 
    geom_sf(aes(fill = lf_mean), colour=NA) +
    theme_void() +
    theme(plot.title = element_text(size = 5)) + 
    annotation_custom(grob = lnd_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                      ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) + 
    scale_fill_distiller(palette = "RdYlBu", limits =c(0, 2), 
                         guide = guide_colourbar(title.position = "top")) + 
    labs(fill=leg_name) + 
    theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.key.height = unit(0.2, 'cm'), #change legend key height
          legend.key.width = unit(bar_width, 'cm'), #change legend key width
          legend.title = element_text(size=7), #change legend title font size
          legend.text = element_text(size=5)) #change legend text font size
}
n_week_plot <- 20
maxy_plot <- 2
rng <- c(0, maxy_plot)
#rev(rev(mid_week_unique)[1 + 1 * (0:(n_week_plot - 1))])
leg_name = "Prevalence"

map_zoom_plot_I <- function(week_plot_curr) {
  IR_week <- IR[IR$mid_week == week_plot_curr, ]
  curr_y <- IR_week[match(LTLA_shp_Reg$LAD20NM, IR_week$ltla), "I_m"]
  curr_y[curr_y > maxy_plot] <- maxy_plot
  curr_y[is.na(curr_y)] <- mean(curr_y, na.rm = T)
  LTLA_shp_Reg$lf_mean <- curr_y
  y_lnd <- ggplot(data = filter(LTLA_shp_Reg, RGN20NM == "London")) + 
    geom_sf(aes(fill = lf_mean ), color = NA, show.legend = FALSE) + 
    theme_void() + ggtitle("London") + 
    scale_fill_viridis_c(option = "A", direction = -1, limits = c(0, 2.5)) + 
    theme(plot.title = element_text(size = 5), 
          legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.key.height = unit(0.2, 'cm'), #change legend key height
          legend.key.width = unit(0.2, 'cm'), #change legend key width
          legend.title = element_text(size=7), #change legend title font size
          legend.text = element_text(size=5)) #change legend text font size
  lnd_grob <- ggplotGrob(y_lnd)
  ggplot(data = LTLA_shp_Reg) + 
    geom_sf(aes(fill = lf_mean), colour=NA) +
    theme_void() +
    ggtitle(week_plot_curr) +
    theme(plot.title = element_text(size = 7)) +
    annotation_custom(grob = lnd_grob, xmin = xmin, xmax = xmin+0.35*(xmax-xmin), 
                      ymin = ymin+3.5*yoffset, ymax =ymin + 3.5*yoffset + 0.35*(ymax-ymin)) + 
    scale_fill_viridis_c(option = "A", direction = -1, limits = c(0, 2.5), 
                         guide = guide_colourbar(title.position = "top")) + 
    labs(fill=leg_name) +
    theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.key.height = unit(0.2, 'cm'), #change legend key height
          legend.key.width = unit(bar_width, 'cm'), #change legend key width
          legend.title = element_text(size=7), #change legend title font size
          legend.text = element_text(size=5)) #change legend text font size
}
ltla_I_plots <- vector("list", 5)
for (i in seq_along(mid_weeks)) {
  mid_week <- mid_weeks[i]
  ltla_I_plots[[i]] <- map_zoom_plot_I(mid_week)
}
IR$week <- ymd(IR$mid_week) - 1

joint_df <- IR %>%
  mutate(mid_week = mid_week + 3) %>%
  filter(mid_week >= as.Date("2021-01-01")) %>%
  select(-week) %>%
  left_join(LTLA2Reg, by = c("ltla" = "LAD20NM")) %>%
  select(-c(ltla)) %>%
  left_join(voc_df, by = c("LAD20CD", "mid_week")) %>%
  #  filter(!is.na(ltla)) %>%
  rename(Region = RGN20NM) %>%
  mutate(Region = factor(Region))
scatt_plots <- vector("list", 5)
for (i in seq_along(mid_weeks)) {
  this_mid_week <- mid_weeks[i] + 3
  scatt_plots[[i]] <- joint_df %>%
    filter(!is.na(R_m)) %>%
    filter(mid_week == this_mid_week) %>%
    ggplot(aes(`Percent delta VOC`, R_m, color = Region)) +
    scale_colour_viridis_d(option = "B", alpha = 0.7, na.translate = FALSE) + 
    geom_point(size = 0.5) + 
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.5) +
    theme_minimal() + 
    labs(y = "Rt") +#+ ggtitle("2020-11-07")+
    theme(plot.title = element_text(size = 5), 
          axis.title = element_text(size = 5),
          axis.text = element_text(size = 5),
          legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.key.height = unit(0.2, 'cm'), #change legend key height
          legend.key.width = unit(0.4, 'cm'), #change legend key width
          legend.title = element_text(size=7), #change legend title font size
          legend.text = element_text(size=5)) #change legend text font size
}
p1 <- ltla_I_plots[[1]]  + ltla_R_plots[[1]] + var_plots[[1]] + scatt_plots[[1]] + theme(aspect.ratio = 1/1) +
  ltla_I_plots[[2]]  + ltla_R_plots[[2]] + var_plots[[2]] + scatt_plots[[2]] + theme(aspect.ratio = 1/1) +
  ltla_I_plots[[3]]  + ltla_R_plots[[3]] + var_plots[[3]] + scatt_plots[[3]] + theme(aspect.ratio = 1/1) +
  ltla_I_plots[[4]]  + ltla_R_plots[[4]] + var_plots[[4]] + scatt_plots[[4]] + theme(aspect.ratio = 1/1) +
  ltla_I_plots[[5]]  + ltla_R_plots[[5]] + var_plots[[5]] + scatt_plots[[5]] + theme(aspect.ratio = 1/1) +
  plot_layout(ncol = 4,
              byrow = TRUE, guides = 'collect') & 
  guides(color = guide_legend(nrow = 9)) &
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=9)) # plot.margin =margin(-1,-1,-1,, "cm"))

plot_dir <- "plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_dir, "epimap_maps.png"), 
       width = 8, height = 12, units = "in")
