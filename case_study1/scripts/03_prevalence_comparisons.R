dir.create("output", showWarnings = FALSE)

####################
# Get REACT dates
####################

react_round_df <- readr::read_csv("data/react_round.csv")
react_rounds <- sort(unique(react_round_df$round))

round_start_dates <- as.Date(c("2020-11-13", "2021-01-06", "2021-02-04",
                               "2021-03-11", "2021-04-15", "2021-05-20",
                               "2021-06-24"))
round_end_dates <- as.Date(c("2020-12-03", "2021-01-22", "2021-02-23",
                             "2021-03-30", "2021-05-03", "2021-06-07",
                             "2021-07-12"))

react_date_df <- tibble(round = react_rounds, 
                        start_date = round_start_dates,
                        end_date = round_end_dates) %>%
  mutate(mid_date = as.Date(round(0.5 * (as.numeric(start_date) + 
                                     as.numeric(end_date))),
                            origin = "1970-01-01")) %>%
  right_join(distinct(ltla_df, mid_week), by = character()) %>%
  group_by(round) %>%
  filter(abs(mid_week - mid_date) == min(abs(mid_week - mid_date)))
readr::write_csv(react_date_df, "data/react_dates.csv")

######################################
# Calculate raw pillar 2 estimates
######################################

raw_pillar2_df <- ltla_df %>%
  inner_join(react_date_df, by = "mid_week") %>%
#  group_by(round) %>%
#  filter(abs(mid_week - mid_date) == min(abs(mid_week - mid_date))) %>%
  bind_cols(as_tibble(Hmisc::binconf(.$nt, .$Nt) * 100)) %>%
  rename(m = PointEst, l = Lower, u = Upper)

readr::write_csv(raw_pillar2_df, "data/raw_pillar2.csv")

############################################
# Calculate REACT LTLA exact Binomial CIs
############################################

old_bucks_names <- c("Aylesbury Vale", "Chiltern", "South Bucks", "Wycombe")
new_bucks_name <- "Buckinghamshire"

old_north_names <- c("Corby", "East Northamptonshire", "Kettering", "Wellingborough")
new_north_name <- "North Northamptonshire"

old_west_names <- c("Daventry", "Northampton", "South Northamptonshire")
new_west_name <- "West Northamptonshire"

react_ltla_df <- react_round_df %>%
  mutate(ltla = case_when(ltla %in% old_bucks_names ~ new_bucks_name,
                          ltla %in% old_north_names ~ new_north_name,
                          ltla %in% old_west_names ~ new_west_name,
                          TRUE ~ ltla)) %>%
  group_by(ltla, round) %>%
  summarise(positive = sum(positive),
            number_samples = sum(number_samples), .groups = "drop") %>%
  bind_cols(as_tibble(Hmisc::binconf(.$positive, .$number_samples) * 100)) %>%
  rename(m = PointEst, l = Lower, u = Upper)

readr::write_csv(react_ltla_df, "data/react_ltla.csv")
