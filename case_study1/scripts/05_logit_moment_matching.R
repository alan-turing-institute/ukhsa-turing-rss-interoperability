ltla_df <- readr::read_csv("data/ltla.csv")
mid_week_unique <- unique(ltla_df$mid_week)

logit_moments <- list()

for (this_ltla in ltla_names) {
  this_M <- ltla_df %>%
    filter(ltla == this_ltla) %>%
    pull(M) %>%
    first()
  
  I_prop <- control$I_seq / this_M
  
  this_mom <- apply(ltla_prevalence[[this_ltla]]$log_lik, 1, logit_mom, I_prop)
  
  logit_moments[[this_ltla]] <- do.call(rbind, this_mom) %>%
    mutate(ltla = this_ltla, mid_week = mid_week_unique, M = this_M)
}

logit_moment_df <- do.call(bind_rows, logit_moments)

readr::write_csv(logit_moment_df, file.path(output_dir, "/logit_moments.csv"))
