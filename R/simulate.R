# Simulating allometric data ----

#' Target covariance matrix for appendage, mass, and an optional gradient
#'
#' Assembles the log-scale covariance matrix implied by an allometry, optionally
#' extended with an environmental gradient. The allometry is resolved by
#' [implied_allometry()], so supply any three of `r_app_mass`, `b_ols`, `b_sma`,
#' `b_avg`, `sd_log_append`, `sd_log_mass`.
#'
#' The gradient enters at **unit standard deviation by default**, because this
#' matrix primarily describes a correlation structure and [sim_allometric()]
#' rescales the drawn gradient itself. Set `sd_gradient` to put the gradient on
#' its own scale — e.g. `sd_gradient = 0.18` for a temperature index — so the
#' returned matrix is the covariance of your actual system rather than a
#' standardised one.
#'
#' @inheritParams implied_allometry
#' @param gradient Optional string naming an environmental gradient (e.g.
#'   `"Temperature"`). `NULL` (the default) returns the 2x2 morphological block.
#' @param r_grad_app,r_grad_mass Correlations of the gradient with
#'   `log(Append)` and `log(Mass)`. Ignored when `gradient` is `NULL`.
#' @param sd_gradient Standard deviation of the gradient. Defaults to `1`.
#'   Ignored when `gradient` is `NULL`. The morphological correlations and
#'   slopes do not depend on it; it scales only the gradient's own variance and
#'   covariances.
#'
#' @return A covariance matrix, 2x2 with dimnames `Append`, `Mass`, or 3x3 with
#'   the gradient's name appended.
#'
#' @seealso [sim_allometric()], which draws from this structure.
#'
#' @examples
#' build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07)
#'
#' # A temperature gradient on its own scale (SD 0.18), not standardised
#' build_cov_mat(
#'   b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07,
#'   gradient = "Temperature", r_grad_app = -0.3, r_grad_mass = -0.1,
#'   sd_gradient = 0.18
#' )
#'
#' @export
build_cov_mat <- function(r_app_mass = NULL, b_ols = NULL, b_sma = NULL, b_avg = NULL,
                          sd_log_append = NULL, sd_log_mass = NULL,
                          gradient = NULL, r_grad_app = 0, r_grad_mass = 0,
                          sd_gradient = 1) {
  parms <- solve_allometry(
    r_app_mass = r_app_mass, b_ols = b_ols, b_sma = b_sma, b_avg = b_avg,
    sd_log_append = sd_log_append, sd_log_mass = sd_log_mass
  )

  if (is.null(gradient)) {
    vars <- c("Append", "Mass")
    cor_mat <- matrix(c(1, parms$r_app_mass, parms$r_app_mass, 1),
                      nrow = 2, dimnames = list(vars, vars))
    check_pos_def(cor_mat)
    return(cor_to_cov(cor_mat, c(parms$sd_log_append, parms$sd_log_mass)))
  }

  check_gradient_name(gradient)
  if (!rlang::is_scalar_double(sd_gradient) && !rlang::is_scalar_integer(sd_gradient) ||
      sd_gradient <= 0) {
    rlang::abort("`sd_gradient` must be a single positive number.")
  }
  vars <- c("Append", "Mass", gradient)
  cor_mat <- matrix(
    c(1,            parms$r_app_mass, r_grad_app,
      parms$r_app_mass, 1,            r_grad_mass,
      r_grad_app,   r_grad_mass,      1),
    nrow = 3, byrow = TRUE, dimnames = list(vars, vars)
  )
  check_pos_def(cor_mat)
  cor_to_cov(cor_mat, c(parms$sd_log_append, parms$sd_log_mass, sd_gradient))
}


### A gradient name becomes a column name, so it must be a single usable string that does not collide with the morphological columns.
check_gradient_name <- function(gradient, call = rlang::caller_env()) {
  if (!rlang::is_string(gradient) || !nzchar(gradient)) {
    rlang::abort("`gradient` must be a single non-empty string naming the gradient column, or `NULL`.", call = call)
  }
  if (gradient %in% c("Append", "Mass", "Append_log", "Mass_log")) {
    rlang::abort(paste0("`gradient` cannot be named `", gradient, "`; that column already exists."), call = call)
  }
  invisible(gradient)
}

### Resolve the gradient's location and spread from either `gradient_range` (its bounds) or an explicit mean and SD, refusing to accept both.
resolve_gradient_scale <- function(gradient_range, mean_gradient, sd_gradient,
                                   gradient_dist, scale_supplied, call) {
  if (is.null(gradient_range)) return(list(mean = mean_gradient, sd = sd_gradient))

  if (scale_supplied) {
    rlang::abort(
      c("Supply either `gradient_range` or `mean_gradient`/`sd_gradient`, not both.",
        i = "`gradient_range` already determines the mean and standard deviation."),
      call = call
    )
  }
  if (gradient_dist != "uniform") {
    rlang::abort(
      c("`gradient_range` requires `gradient_dist = \"uniform\"`.",
        i = "A normal gradient is unbounded, so it has no range."),
      call = call
    )
  }
  if (!is.numeric(gradient_range) || length(gradient_range) != 2L ||
      anyNA(gradient_range) || gradient_range[1] >= gradient_range[2]) {
    rlang::abort("`gradient_range` must be two increasing numbers, e.g. `c(1970, 2026)`.", call = call)
  }
  list(
    mean = mean(gradient_range),
    sd   = diff(gradient_range) / sqrt(12)   # SD of Uniform(a, b)
  )
}

### Draw the standardised gradient. Uniform on +/-sqrt(3) has mean 0 and SD 1, matching the standard normal, so the two marginals are interchangeable in the construction below.
# Exact correlations require the drawn gradient to have sample variance exactly 1, which forces a z-score. Z-scoring an iid uniform draw inflates its extremes past +/-sqrt(3), so `gradient_range = c(-90, 90)` would emit latitudes beyond +/-90. A systematic sample is exactly symmetric, so its z-score stays strictly inside sqrt(3)*sqrt((n-1)/(n+1)) < sqrt(3), and bounds hold. Order is shuffled so no row position carries gradient information.
draw_gradient <- function(n, gradient_dist, empirical) {
  if (!empirical) {
    return(switch(
      gradient_dist,
      uniform = stats::runif(n, -sqrt(3), sqrt(3)),
      normal  = stats::rnorm(n)
    ))
  }
  if (gradient_dist == "uniform") {
    g <- stats::qunif(stats::ppoints(n), -sqrt(3), sqrt(3))
    return(sample((g - mean(g)) / stats::sd(g)))
  }
  g <- stats::rnorm(n)
  (g - mean(g)) / stats::sd(g)
}

### Draw residuals with an exact target covariance that are also exactly orthogonal to the gradient. Residualising against `g` removes any chance correlation; whitening then re-imposes the target covariance on the residualised columns, which remain orthogonal to `g` because each is a linear combination of columns that already were.
draw_residuals <- function(n, rho_e, g, empirical) {
  sigma_e <- matrix(c(1, rho_e, rho_e, 1), nrow = 2)
  e <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma_e, empirical = empirical)
  if (!empirical) return(e)

  e <- stats::residuals(stats::lm(e ~ g))
  e %*% solve(chol(stats::cov(e))) %*% chol(sigma_e)
}


#' Simulate log-normal morphology, optionally along an environmental gradient
#'
#' Draws appendage length and body mass from a bivariate log-normal with a
#' specified allometry, and optionally an exogenous environmental gradient
#' (temperature, latitude, year, rainfall) that both traits respond to linearly.
#'
#' The allometry is resolved by [implied_allometry()]. Supply **two** of
#' `r_app_mass`, `b_ols`, `b_sma`, `b_avg`, or none to accept the isometric
#' default `b_sma = 1/3, r_app_mass = 0.3`. Supplying exactly one is an error.
#' `sd_log_append` and `sd_log_mass` set the scale and default to
#' `sd_log_mass = 0.07`; they leave every slope and correlation unchanged.
#'
#' @section The gradient:
#' The gradient is a **sampling design variable**, not a third trait, so it is
#' drawn exogenously and the log traits are built linearly from it. Its marginal
#' is therefore free, and it defaults to `"uniform"` because:
#'
#' - **it is bounded**, so `gradient_range = c(-90, 90)` really does yield
#'   latitudes within \eqn{\pm 90}, and years stay inside their range;
#' - **it covers the gradient evenly.** Across 15 equal-width bins a normal
#'   gradient leaves the extreme bins nearly empty, with the fullest bin
#'   thousands of times the emptiest; a uniform one is flat.
#'
#' `gradient_dist = "normal"` recovers the ordinary three-variable multivariate
#' normal. Either way `r_grad_app` is the correlation of the gradient with
#' `log(Append)` (not with `Append`), and \eqn{E[\log A \mid G]} is linear in
#' \eqn{G} by construction.
#'
#' Not every set of correlations is attainable. An appendage that lengthens
#' while mass falls along the gradient is incompatible with a strong positive
#' appendage-mass correlation, and such a request is rejected rather than
#' silently approximated.
#'
#' @section Error:
#' Error is added on the raw scale, after exponentiation, each component a
#' fraction of the trait's own standard deviation. `meas_error` is observer
#' error and applies to both traits alike; `transient_error_*` represents real
#' short-term biological fluctuation, which afflicts mass far more than a
#' skeletal appendage.
#'
#' @section Reproducibility and trimming:
#' This function consumes the random number stream, so call [set.seed()] before
#' it for reproducible output. Because the internal order of the draws is an
#' implementation detail, byte-for-byte reproducibility is tied to a fixed
#' package version: pin one (e.g. `remotes::install_github("LezzGitIt/sliR@v0.1.0")`)
#' for analyses that must reproduce exactly.
#'
#' Trimming (`trim_sd`) is applied **jointly** to `Append` and `Mass`: a row is
#' dropped if either trait exceeds the cut. Code that instead trims one trait
#' and then the other, using each trait's pre-trim standard deviation, will
#' retain a slightly different set of rows for the same draw.
#'
#' @inheritParams build_cov_mat
#' @param n Number of individuals to draw, before trimming.
#' @param mean_append,mean_mass Raw-scale means, used as the log-scale location
#'   via `log()`.
#' @param gradient_dist Marginal distribution of the gradient: `"uniform"` (the
#'   default) or `"normal"`.
#' @param gradient_range Two increasing numbers giving the gradient's bounds,
#'   e.g. `c(1970, 2026)`. Uniform gradients only. Determines `mean_gradient`
#'   and `sd_gradient`, which must then not be supplied.
#' @param mean_gradient,sd_gradient Location and spread of the gradient.
#'   Default to `0` and `1`, a z-scored gradient.
#' @param meas_error Observer error, as a fraction of each trait's standard
#'   deviation. Applied to both traits.
#' @param transient_error_append,transient_error_mass Biological fluctuation, as
#'   a fraction of that trait's standard deviation.
#' @param trim_sd Drop individuals more than `trim_sd` raw-scale standard
#'   deviations from the mean of `Append` or `Mass`. `NULL` disables trimming.
#'   The gradient is never trimmed. Because the traits are log-normal, raw-scale
#'   trimming is asymmetric and removes more of the right tail than the left.
#' @param log If `TRUE` (the default), append the `Append_log` and `Mass_log`
#'   columns, computed after error is added. `FALSE` returns the raw traits and
#'   the gradient only.
#' @param empirical If `TRUE` (the default), the drawn sample has exactly the
#'   requested moments and correlations, rather than being a random draw from a
#'   population with them. Trimming and error break this guarantee.
#'
#' @return A tibble with `n` rows or fewer: `Append`, `Mass`, the gradient
#'   column if requested, and `Append_log`, `Mass_log` unless `log = FALSE`.
#'
#' @seealso [implied_allometry()] to inspect what a parameter set implies, and
#'   [sim_correlated()] for a raw-scale bivariate draw with no allometric target.
#'
#' @examples
#' set.seed(1)
#'
#' # Isometry, no gradient
#' d <- sim_allometric(n = 500)
#' coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]]
#'
#' # Allen's rule: warmer places, relatively longer appendages
#' warm <- sim_allometric(
#'   n = 500, gradient = "Temperature", gradient_range = c(0, 24),
#'   r_grad_app = 0.3, r_grad_mass = -0.1
#' )
#' range(warm$Temperature)
#'
#' # A century of specimens, drawn evenly across years
#' yr <- sim_allometric(n = 500, gradient = "Year", gradient_range = c(1970, 2026))
#'
#' # Raw traits only
#' names(sim_allometric(n = 10, log = FALSE))
#'
#' @export
sim_allometric <- function(n = 3000,
                           r_app_mass = NULL, b_ols = NULL, b_sma = NULL, b_avg = NULL,
                           sd_log_append = NULL, sd_log_mass = NULL,
                           mean_append = 180, mean_mass = 80,
                           gradient = NULL,
                           gradient_dist = c("uniform", "normal"),
                           gradient_range = NULL,
                           mean_gradient = 0, sd_gradient = 1,
                           r_grad_app = 0, r_grad_mass = 0,
                           meas_error = 0,
                           transient_error_append = 0, transient_error_mass = 0,
                           trim_sd = 3, log = TRUE, empirical = TRUE) {
  ## `missing()` must be consulted before match.arg() assigns to `gradient_dist`, or it always reports FALSE.
  scale_supplied <- !missing(mean_gradient) || !missing(sd_gradient)
  grad_supplied  <- scale_supplied || !is.null(gradient_range) ||
    !missing(r_grad_app) || !missing(r_grad_mass) || !missing(gradient_dist)
  gradient_dist  <- match.arg(gradient_dist)

  parms <- solve_allometry(
    r_app_mass = r_app_mass, b_ols = b_ols, b_sma = b_sma, b_avg = b_avg,
    sd_log_append = sd_log_append, sd_log_mass = sd_log_mass
  )
  mu_a <- base::log(mean_append)
  mu_m <- base::log(mean_mass)

  if (is.null(gradient)) {
    if (grad_supplied) {
      rlang::abort(
        c("Gradient arguments were supplied without naming a gradient.",
          i = "Set `gradient = \"Temperature\"` (or another name) to include one."))
    }
    Sigma <- cor_to_cov(
      matrix(c(1, parms$r_app_mass, parms$r_app_mass, 1), nrow = 2),
      c(parms$sd_log_append, parms$sd_log_mass)
    )
    draw <- MASS::mvrnorm(n, mu = c(mu_a, mu_m), Sigma = Sigma, empirical = empirical)
    sim  <- tibble::tibble(Append = exp(draw[, 1]), Mass = exp(draw[, 2]))
  } else {
    check_gradient_name(gradient)
    ## Reject impossible correlation triples up front. This is the same determinant condition that makes the residual covariance below positive definite, so one check serves both.
    vars <- c("Append", "Mass", gradient)
    check_pos_def(matrix(
      c(1,                parms$r_app_mass, r_grad_app,
        parms$r_app_mass, 1,                r_grad_mass,
        r_grad_app,       r_grad_mass,      1),
      nrow = 3, byrow = TRUE, dimnames = list(vars, vars)
    ))

    g_scale <- resolve_gradient_scale(gradient_range, mean_gradient, sd_gradient,
                                      gradient_dist, scale_supplied, rlang::current_env())

    g <- draw_gradient(n, gradient_dist, empirical)

    ## Residual correlation chosen so that cor(log_A, log_M) lands on r_app_mass once the shared gradient path is added back in.
    rho_e <- (parms$r_app_mass - r_grad_app * r_grad_mass) /
      sqrt((1 - r_grad_app^2) * (1 - r_grad_mass^2))
    e <- draw_residuals(n, rho_e, g, empirical)

    log_a <- mu_a + parms$sd_log_append * (r_grad_app  * g + sqrt(1 - r_grad_app^2)  * e[, 1])
    log_m <- mu_m + parms$sd_log_mass   * (r_grad_mass * g + sqrt(1 - r_grad_mass^2) * e[, 2])

    sim <- tibble::tibble(Append = exp(log_a), Mass = exp(log_m))
    sim[[gradient]] <- g_scale$mean + g_scale$sd * g
  }

  sim <- trim_outliers(sim, cols = c("Append", "Mass"), n_sd = trim_sd)
  sim <- dplyr::mutate(
    sim,
    Append = add_error(Append, meas_error, transient_error_append),
    Mass   = add_error(Mass,   meas_error, transient_error_mass)
  )
  if (log) {
    sim <- dplyr::mutate(sim, Append_log = base::log(Append), Mass_log = base::log(Mass))
  }
  sim
}


#' Simulate [sim_allometric()] across a grid of parameters
#'
#' A thin wrapper for the common case of comparing many scenarios: each row of
#' `params` becomes one call to [sim_allometric()], its columns matched to
#' `sim_allometric()`'s arguments by name exactly as `purrr::pmap()` would.
#' Arguments that are constant across every scenario — most often `gradient`,
#' a single string naming one column — are passed via `...` instead of
#' appearing in `params`.
#'
#' @param params A data frame, one row per scenario, with columns named after
#'   [sim_allometric()]'s arguments (e.g. `n`, `r_app_mass`, `r_grad_app`).
#' @param ... Arguments passed to every call of [sim_allometric()], constant
#'   across scenarios (e.g. `gradient = "Temperature"`).
#' @param id_col Name of the integer column identifying which row of `params`
#'   produced each simulated row, matching that row's position in `params`.
#'   Defaults to `".scenario"`.
#'
#' @return A single tibble: the row-bound output of [sim_allometric()] for
#'   every row of `params`, with `id_col` added. Join `params` back on by row
#'   number (`dplyr::mutate(params, .scenario = dplyr::row_number())`) to
#'   recover the parameter values behind each scenario.
#'
#' @seealso [sim_allometric()], which this maps over.
#'
#' @examples
#' grid <- expand.grid(n = 200, r_app_mass = c(0.1, 0.3, 0.5), b_sma = 1 / 3)
#' set.seed(1)
#' sims <- sim_grid(grid)
#' dplyr::count(sims, .scenario)
#'
#' # Arguments constant across the grid, like the gradient's name and range,
#' # are passed once via `...` rather than repeated in every row of `params`.
#' grad_grid <- expand.grid(n = 200, r_grad_app = c(-0.3, 0, 0.3), r_grad_mass = -0.1)
#' set.seed(1)
#' grad_sims <- sim_grid(grad_grid, gradient = "Temperature", gradient_range = c(0, 24))
#'
#' @export
sim_grid <- function(params, ..., id_col = ".scenario") {
  if (!is.data.frame(params) || nrow(params) == 0L) {
    rlang::abort("`params` must be a data frame with at least one row.")
  }
  if (!rlang::is_string(id_col) || !nzchar(id_col)) {
    rlang::abort("`id_col` must be a single non-empty string.")
  }
  if (id_col %in% names(params)) {
    rlang::abort(paste0("`id_col` (\"", id_col, "\") already names a column of `params`; choose another."))
  }

  sims <- purrr::pmap(params, sim_allometric, ...)
  dplyr::bind_rows(purrr::imap(sims, \(sim, i) dplyr::mutate(sim, "{id_col}" := as.integer(i), .before = 1)))
}


#' Simulate correlated appendage and mass on the raw scale
#'
#' Draws a bivariate normal pair of appendage length and body mass with a given
#' correlation, then optionally perturbs each with measurement error and
#' transient biological fluctuation. Unlike [sim_allometric()] there is no log
#' transform, no allometric target, and no gradient: this is the tool for
#' building intuition about how a fitted slope responds to correlation and to
#' error added on one trait but not the other.
#'
#' @param n Number of individuals.
#' @param r Correlation between appendage length and mass.
#' @param mu_append,mu_mass Means of the two traits.
#' @param sd_append,sd_mass Standard deviations of the two traits.
#' @inheritParams sim_allometric
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
