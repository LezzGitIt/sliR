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

### Four well-populated Age x Sex cells (n = 30 each, r_app_mass = 0.6 -- exact,
### via sim_allometric()'s empirical = TRUE default -- so each cell clears
### n_min_cell/n_min_marginal = 10/15 on its own, each with its own distinct
### b_sma so a per-cell fit is verifiably not just falling through to a
### coarser level) plus one deliberately orphaned cell (Age = "Rare",
### Sex = "X", n = 4): its own cell, its Age marginal, and its Sex marginal
### are all exactly that same n = 4, so all three fail n_min_cell/n_min_marginal
### identically and it must fall all the way to the pooled slope. Every cell
### shares the same r_app_mass, so a pass/fail split is driven only by n
### against n_min_cell = 10 / n_min_marginal = 15, never by r or p -- the
### tests below always pass those two thresholds explicitly rather than
### relying on the (larger) package defaults.
make_cascade_fixture <- function(seed = 11) {
  set.seed(seed)
  cells <- list(
    list(Age = "Adult", Sex = "F", n = 30, b_sma = 0.30),
    list(Age = "Adult", Sex = "M", n = 30, b_sma = 0.40),
    list(Age = "Juv",   Sex = "F", n = 30, b_sma = 0.50),
    list(Age = "Juv",   Sex = "M", n = 30, b_sma = 0.60),
    list(Age = "Rare",  Sex = "X", n = 4,  b_sma = 0.90)
  )
  purrr::map(cells, \(cl) {
    d <- sim_allometric(n = cl$n, b_sma = cl$b_sma, r_app_mass = 0.6, sd_log_mass = 0.07, trim_sd = NULL)
    d$Age <- cl$Age
    d$Sex <- cl$Sex
    d
  }) |>
    purrr::list_rbind()
}

### Same shape, but r_app_mass near zero everywhere, so even the pooled fit fails reliability.
make_unreliable_fixture <- function(n = 200, seed = 11) {
  set.seed(seed)
  d <- sim_allometric(n = n, b_sma = 1/3, r_app_mass = 0.02, sd_log_mass = 0.07, trim_sd = NULL)
  d$Sex <- rep(c("F", "M"), length.out = nrow(d))
  d
}
