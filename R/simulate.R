# Simulating allometric data ----

#' Log-scale covariance matrix for appendage, mass, and temperature
#'
#' Builds the 3 x 3 covariance matrix consumed by [sim_allometric()]. The
#' matrix is parameterised by a target allometric slope rather than by the two
#' morphological standard deviations directly, because in practice one knows
#' the scaling exponent and one trait's variability, not both variabilities.
#'
#' `b_avg` is the mean of the OLS and SMA slopes of `log(Append)` on
#' `log(Mass)`. Since \eqn{b_{SMA} = b_{OLS} / r}, fixing `b_avg` and `r_am`
#' pins the OLS slope at
#' \deqn{b_{OLS} = \frac{2\, b_{avg}\, r_{am}}{r_{am} + 1}}
#' and hence the ratio of the two morphological standard deviations, since
#' \eqn{b_{OLS} = r_{am}\, \sigma_A / \sigma_M}. `vary` chooses which of the
#' two you supply through `sd_log_morph`; the other is solved for.
#'
#' @param b_avg Target average of the OLS and SMA slopes of log appendage on
#'   log mass. The default `0.33` corresponds to geometric similarity.
#' @param r_am,r_at,r_mt Pearson correlations on the log scale between
#'   appendage and mass, appendage and temperature, and mass and temperature.
#'   `r_am` must be non-zero.
#' @param sd_log_morph Standard deviation on the log scale of whichever
#'   morphological trait `vary` names. To reproduce a randomly drawn value, as
#'   opposed to a fixed one, pass e.g. `sd_log_morph = runif(1, 0.05, 0.09)`.
#' @param vary Which trait `sd_log_morph` refers to: `"mass"` (the default) or
#'   `"append"`. The other trait's standard deviation is derived from `b_avg`
#'   and `r_am`.
#' @param sd_temp Standard deviation of temperature.
#'
#' @return A 3 x 3 covariance matrix with dimnames `Append`, `Mass`, `Temp`.
#'
#' @examples
#' build_cov_mat(b_avg = 0.33, r_am = 0.5, r_at = -0.1, r_mt = -0.4)
#'
#' @export
build_cov_mat <- function(b_avg = 0.33,
                          r_am = 0.3,
                          r_at = -0.1,
                          r_mt = -0.1,
                          sd_log_morph = 0.07,
                          vary = c("mass", "append"),
                          sd_temp = 0.18) {
  vary <- match.arg(vary)
  if (isTRUE(all.equal(r_am, 0))) {
    rlang::abort("`r_am` must be non-zero: with no mass-appendage correlation the allometric slope is undefined.")
  }

  b_ols <- (2 * b_avg * r_am) / (r_am + 1)
  if (vary == "append" && isTRUE(all.equal(b_ols, 0))) {
    rlang::abort("`b_avg` must be non-zero when `vary = \"append\"`, otherwise mass has no implied variance.")
  }

  if (vary == "mass") {
    sd_log_mass   <- sd_log_morph
    sd_log_append <- abs(b_ols / r_am * sd_log_mass)
  } else {
    sd_log_append <- sd_log_morph
    sd_log_mass   <- abs(r_am / b_ols * sd_log_append)
  }

  vars   <- c("Append", "Mass", "Temp")
  cor_mat <- matrix(
    c(1,    r_am, r_at,
      r_am, 1,    r_mt,
      r_at, r_mt, 1),
    nrow = 3, byrow = TRUE, dimnames = list(vars, vars)
  )
  check_pos_def(cor_mat)

  cor_to_cov(cor_mat, c(sd_log_append, sd_log_mass, sd_temp))
}


#' Simulate log-normal morphological data with a target allometric slope
#'
#' Draws appendage length, body mass, and temperature from a multivariate
#' normal on the log scale (temperature stays on its natural scale), then
#' exponentiates the morphological traits. The covariance structure comes from
#' [build_cov_mat()], so the realised allometric slope matches `b_avg` up to
#' the distortion introduced by trimming and by any error you add.
#'
#' Error is added on the raw scale, after exponentiation, and each component is
#' expressed as a fraction of the trait's own standard deviation. The two
#' components are statistically identical draws with different interpretations:
#' `meas_error` is observer error, applied to both traits alike, while
#' `transient_error_*` represents real short-term biological fluctuation, which
#' typically afflicts mass far more than a skeletal appendage.
#'
#' @param n Number of individuals to draw, before trimming.
#' @inheritParams build_cov_mat
#' @param mean_mass,mean_append Means of mass and appendage length on the raw
#'   scale; used as the log-scale location via `log()`.
#' @param mean_temp Mean temperature, on its natural scale.
#' @param meas_error Measurement error, as a fraction of each trait's standard
#'   deviation. Applied to both traits.
#' @param transient_error_mass,transient_error_append Biological fluctuation,
#'   as a fraction of that trait's standard deviation.
#' @param trim_sd Drop individuals lying more than `trim_sd` raw-scale standard
#'   deviations from the mean of `Append` or `Mass`. `NULL` disables trimming.
#'   Because the traits are log-normal, raw-scale trimming is asymmetric and
#'   removes more of the right tail than the left.
#' @param empirical If `TRUE` (the default), the drawn sample has exactly the
#'   specified means and covariance, rather than being a random draw from a
#'   population with those parameters. Trimming and error break this guarantee.
#'
#' @return A tibble with `n` rows or fewer: `Append`, `Mass`, `Temp`, and the
#'   log-scale `Append_log` and `Mass_log`, computed after error is added.
#'
#' @seealso [sim_correlated()] for a simpler raw-scale bivariate draw with no
#'   allometric target.
#'
#' @examples
#' set.seed(1)
#' d <- sim_allometric(n = 500, b_avg = 0.33, r_am = 0.3)
#' coef(smatr::sma(Append_log ~ Mass_log, data = d))
#'
#' # Mass fluctuates transiently; the appendage does not
#' set.seed(1)
#' sim_allometric(n = 500, transient_error_mass = 0.5)
#'
#' @export
sim_allometric <- function(n = 3000,
                           b_avg = 0.33,
                           r_am = 0.3,
                           r_at = -0.1,
                           r_mt = -0.1,
                           mean_mass = 80,
                           mean_append = 180,
                           mean_temp = 1,
                           sd_log_morph = 0.07,
                           vary = c("mass", "append"),
                           sd_temp = 0.18,
                           meas_error = 0,
                           transient_error_mass = 0,
                           transient_error_append = 0,
                           trim_sd = 3,
                           empirical = TRUE) {
  vary  <- match.arg(vary)
  Sigma <- build_cov_mat(
    b_avg = b_avg, r_am = r_am, r_at = r_at, r_mt = r_mt,
    sd_log_morph = sd_log_morph, vary = vary, sd_temp = sd_temp
  )
  mu <- c(log(mean_append), log(mean_mass), mean_temp)

  sim_log <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma, empirical = empirical)
  colnames(sim_log) <- c("Append", "Mass", "Temp")

  sim <- tibble::tibble(
    Append = exp(sim_log[, "Append"]),
    Mass   = exp(sim_log[, "Mass"]),
    Temp   = sim_log[, "Temp"]
  )
  sim <- trim_outliers(sim, cols = c("Append", "Mass"), n_sd = trim_sd)

  sim |>
    dplyr::mutate(
      Append = add_error(Append, meas_error, transient_error_append),
      Mass   = add_error(Mass,   meas_error, transient_error_mass),
      Append_log = log(Append),
      Mass_log   = log(Mass)
    )
}


#' Simulate correlated appendage and mass on the raw scale
#'
#' Draws a bivariate normal pair of appendage length and body mass with a given
#' correlation, then optionally perturbs each with measurement error and
#' transient biological fluctuation. Unlike [sim_allometric()] there is no log
#' transform and no allometric target: this is the tool for building intuition
#' about how a fitted slope responds to correlation and to error added on one
#' trait but not the other.
#'
#' @param n Number of individuals.
#' @param r Correlation between appendage length and mass.
#' @param mu_append,mu_mass Means of the two traits.
#' @param sd_append,sd_mass Standard deviations of the two traits.
#' @inheritParams sim_allometric
#' @param empirical If `TRUE` (the default), the sample has exactly the
#'   specified moments before error is added.
#'
#' @return A tibble of `n` rows with columns `Append` and `Mass`.
#'
#' @examples
#' set.seed(1)
#' # Transient fluctuation in mass alone attenuates the fitted slope
#' clean <- sim_correlated(r = 0.3, transient_error_mass = 0)
#' noisy <- sim_correlated(r = 0.3, transient_error_mass = 1)
#' c(clean = cor(clean$Append, clean$Mass), noisy = cor(noisy$Append, noisy$Mass))
#'
#' @export
sim_correlated <- function(n = 3000,
                           r = 0.3,
                           mu_append = 180,
                           mu_mass = 80,
                           sd_append = 10,
                           sd_mass = 5,
                           meas_error = 0,
                           transient_error_append = 0,
                           transient_error_mass = 0,
                           empirical = TRUE) {
  vars    <- c("Append", "Mass")
  cor_mat <- matrix(c(1, r, r, 1), nrow = 2, dimnames = list(vars, vars))
  check_pos_def(cor_mat)
  Sigma <- cor_to_cov(cor_mat, c(sd_append, sd_mass))

  sim <- MASS::mvrnorm(n, mu = c(mu_append, mu_mass), Sigma = Sigma, empirical = empirical)
  colnames(sim) <- vars

  tibble::as_tibble(sim) |>
    dplyr::mutate(
      Append = add_error(Append, meas_error, transient_error_append),
      Mass   = add_error(Mass,   meas_error, transient_error_mass)
    )
}
