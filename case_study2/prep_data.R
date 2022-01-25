library(readr)
library(dplyr)
library(readxl)
library(tidyr)
library(readODS)
library(lubridate)

## space
library(rgdal)
library(sf)
library(sp)
library(raster)
library(spdep)

start_date <- as.Date("2020-05-31")
end_date <- as.Date("2021-06-20")

# Islands (Wight and Scilly)
island_codes <- c("E06000046", "E06000053")

merge_Hackney_CoL <- TRUE
#filter_REG <- NULL

# 1. Space data -----------------------------------------------------------

## Local Authority Districts (same as LTLA, checked) to Regions lookup:
LTLA2Reg <- read_csv("data/geography/Local_Authority_District_to_Region_(April_2021)_Lookup_in_England.csv")

# import LTLA shapefile
LTLA_shp = readOGR("data/geography/Local_Authority_Districts_(May_2021)_UK_BFE_V3/Local_Authority_Districts_(May_2021)_UK_BFE_V3.shp")
LTLA_shp <- st_as_sf(LTLA_shp)

# remove LTLAs for which we don't have testing data (i.e. all that are not in England)
LTLA_shp_EN = LTLA_shp %>% dplyr::filter(LAD21CD %in% LTLA2Reg$LAD21CD)

# remove island
islands = c("Isle of Wight",   "Isles of Scilly")
LTLA_shp_EN = LTLA_shp_EN[!LTLA_shp_EN$LAD21NM %in% islands,]

crs <- st_crs(LTLA_shp_EN)

# merge together Hackney and City of London 
if (merge_Hackney_CoL) {
  
  # assign to City of London the code of Hackney and use it as
  # grouping variable
  LTLA_shp_EN$LAD21CD[LTLA_shp_EN$LAD21NM == "City of London"] <- "E09000012"
  
  st_union_by = function(geo, group) {
    
    y2 = list()
    #loop over by groups and merge units
    for (i in unique(group)) {
      #which units
      z = geo[group == i]
      #merge
      y = Reduce(st_union, z)
      y2[[i]] = y
    }
    
    y3 <- st_sfc(y2)
    #
    st_sf(data.frame(LAD21CD = names(y3), geom=st_sfc(y3)))
  }
  
  LTLA_shp_EN <- st_union_by(LTLA_shp_EN$geometry, LTLA_shp_EN$LAD21C)
}

# add region lookup to LTLA shapefile
LTLA_shp_Reg  <- left_join(LTLA_shp_EN, LTLA2Reg, by = c("LAD21CD" = "LAD21CD"))

# add area index 
LTLA_shp_Reg$LTLA_ID <- 1:nrow(LTLA_shp_Reg)

# Determine neighborhood/adjacency information needed for neighborhood-based CAR model
W.nb <- spdep::poly2nb(LTLA_shp_Reg)
nb2INLA("data/W.adj", W.nb)

st_crs(LTLA_shp_Reg) <- crs
save(LTLA_shp_Reg, file="data/space_obj.RData")

# Prevalence --------------------------------------------------------------

year_week <- function(x,base) isoweek(x) - isoweek(base) + 53*(year(x) - year(base))

raw_prev_moments <- read_csv("data/logit_moments.csv", col_types = cols(mid_week = col_date(format = "%Y-%m-%d")))

# Get index for each week/year
raw_prev_moments$week <- year_week(raw_prev_moments$mid_week, start_date)
# fix end of year week
raw_prev_moments <- raw_prev_moments %>% mutate(week = if_else(week>80, week-53, week))

prev_df <- raw_prev_moments %>% 
  filter(between(mid_week, start_date, end_date)) %>%
  mutate(AppointmentDate = week)

prev_data <- LTLA_shp_Reg %>%
  left_join(prev_df, by = c("LAD21NM" = "ltla")) %>%
  mutate(LAD20NM = factor(LAD21NM, levels = unique(LAD21NM))) 

prev_data$geometry <- NULL

prev_data <- prev_data %>% 
  mutate(date_ID = as.integer(as.factor(AppointmentDate)))

set.seed(111)
prev_data <- prev_data %>%  # left_join(prev_data, n_pop_mat) %>%
  mutate(counts = round(boot::inv.logit(mean)*M) + rpois(n = nrow(prev_data), lambda = 1)) 

#sum(is.na(prev_data))


# BAME proportion ---------------------------------------------------------

bame_raw <- read_csv("data/census2011_ethnicity_LSOA_black_southasian.csv", skip = 7)

#LSOA2LTLA <-  read_csv("data/geography/PCD_OA_LSOA_MSOA_LAD_AUG21_UK_LU.csv", 
#                       col_types = cols(pcd7 = col_skip(), pcd8 = col_skip(), 
#                                        pcds = col_skip(), dointr = col_skip(), 
#                                        doterm = col_skip(), usertype = col_skip(), 
#                                        oa11cd = col_skip())) %>% distinct()
#LSOA2LTLA <- distinct(LSOA2LTLA)
#write.csv(LSOA2LTLA, 'data/geography/LSOA2LTLA.csv')

LSOA2LTLA <- read_csv("data/geography/LSOA2LTLA.csv")

bame_raw <- bame_raw %>% 
  rename(LSOA11CD = `2011 super output area - lower layer`) %>% 
  mutate(LSOA11CD = gsub(" :.*","", LSOA11CD))

bame_prop <-  LSOA2LTLA %>% 
  left_join(bame_raw, by = c( "lsoa11cd" = "LSOA11CD" )) %>% 
  group_by(LAD21CD = ladcd) %>% mutate(Black = `Black/African/Caribbean/Black British`, 
                               South_Asian = `Asian/Asian British: Pakistani` + `Asian/Asian British: Indian` + `Asian/Asian British: Bangladeshi`,
                               Other_BAME = `All usual residents` - White - `Black/African/Caribbean/Black British` - 
                                 (`Asian/Asian British: Pakistani` + `Asian/Asian British: Indian` + 
                                    `Asian/Asian British: Bangladeshi`)) %>%
  summarise_at(c("All usual residents", "Black", "South_Asian", "White", "Other_BAME"), sum) %>%
  mutate(Black_prop = Black/`All usual residents`, 
         South_Asian_prop = South_Asian/`All usual residents`,
         Other_BAME_prop = Other_BAME/`All usual residents`,
         BAME = 1-White/`All usual residents`) %>% 
  dplyr::select(LAD21CD, Black_prop,South_Asian_prop,Other_BAME_prop, BAME )

bame_prop <- bame_prop %>% ungroup %>% mutate(lad21cd = LAD21CD, 
                                              bame_quint = gtools::quantcut(BAME, q=5, na.rm=TRUE))


prev_data = left_join(prev_data, bame_prop)

# IMD ---------------------------------------------------------------------

# https://www.gov.uk/government/statistics/english-indices-of-deprivation-2019

IMD_LTLA_2019 <- read_excel("data/File_10_-_IoD2019_Local_Authority_District_Summaries__lower-tier__.xlsx", 
                            sheet = "IMD", n_max = 318) %>% rename(lad21cd = `Local Authority District code (2019)`) %>%
                              mutate(lad21cd = if_else(lad21cd %in%c("E07000004", "E07000005", "E07000006", "E07000007"),
                                                                                     "E06000060", lad21cd)) %>%
                             mutate(lad21cd = if_else(lad21cd %in%c("E07000150", "E07000152", "E07000153", "E07000156"),
                                                      "E06000061", lad21cd)) %>%
                            mutate(lad21cd = if_else(lad21cd %in%c("E07000151", "E07000154", "E07000155"),
                                                      "E06000062", lad21cd)) %>%
  group_by(lad21cd) %>% summarise(IMD = mean(`IMD - Average score`)) %>%
  mutate(IMD_quint = factor(ntile(IMD, 5)))


prev_data <- left_join(prev_data, IMD_LTLA_2019, by = c("lad21cd" = "lad21cd"))
prev_data$IMD_stand <- scale(prev_data$IMD)
prev_data$bame_stand <- scale(prev_data$BAME)

saveRDS(prev_data, "prev_data.rds")

