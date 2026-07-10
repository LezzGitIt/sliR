### Two age classes with deliberately different allometric slopes, and a sex label assigned independently of age. The Sex model therefore recovers a slope near the pooled value while the Age model splits, so every age x sex cell averages two slopes that disagree — the situation `slope_diff_warn` exists to catch.
make_divergent_slopes <- function(n_per_age = 300, seed = 11) {
  set.seed(seed)
  shallow <- sim_allometric(n = n_per_age, b_sma = 0.25, r_app_mass = 0.3, sd_log_mass = 0.07, trim_sd = NULL)
  steep   <- sim_allometric(n = n_per_age, b_sma = 0.85, r_app_mass = 0.3, sd_log_mass = 0.07, trim_sd = NULL)
  shallow$Age <- "Juv"
  steep$Age   <- "Adult"
  d <- rbind(shallow, steep)
  d$Sex <- rep(c("F", "M"), length.out = nrow(d))
  d
}

### Same structure, but both age classes share one slope, so the Age and Sex models agree and no warning should fire.
make_concordant_slopes <- function(n_per_age = 300, seed = 11) {
  set.seed(seed)
  d <- sim_allometric(n = 2 * n_per_age, b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07, trim_sd = NULL)
  d$Age <- rep(c("Juv", "Adult"), each = 2, length.out = nrow(d))
  d$Sex <- rep(c("F", "M"), length.out = nrow(d))
  d
}
