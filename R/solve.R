# Allometry solver ----

### A log-log allometry has three free parameters: the two log-scale standard deviations and the correlation between them. Everything else is derived, via b_sma = sd_A/sd_M (signed by r), b_ols = r * b_sma, and b_avg = (b_ols + b_sma)/2. So any two of the four shape quantities determine the other two, and one standard deviation then fixes the scale.
shape_args <- c("r_app_mass", "b_ols", "b_sma", "b_avg")
scale_args <- c("sd_log_append", "sd_log_mass")

### The trio applied when the user supplies no allometry at all. b_sma = 1/3 is isometry as an SMA fit would report it, which is the exponent calc_sli() assumes.
default_allometry <- list(b_sma = 1 / 3, r_app_mass = 0.3, sd_log_mass = 0.07)

check_scalar_number <- function(x, nm, call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    rlang::abort(paste0("`", nm, "` must be a single non-missing number."), call = call)
  }
  invisible(x)
}

abort_under_determined <- function(shape, call) {
  rlang::abort(
    c(
      paste0("The allometry is under-determined: only `", shape,
             "` was supplied, which fixes one of its two shape quantities."),
      i = paste0("Add a second, from: ",
                 paste0("`", setdiff(shape_args, shape), "`", collapse = ", "),
                 ". Or supply both standard deviations."),
      x = "Pairing your value with a default correlation would silently determine the remaining slopes.",
      i = "Supply no shape arguments at all to accept the isometric default (`b_sma = 1/3, r_app_mass = 0.3`)."
    ),
    call = call
  )
}

### Resolve the correlation and the SMA slope, from which every other shape quantity follows.
solve_shape <- function(s, shape, scale, call) {
  nonzero <- function(x, nm) {
    if (isTRUE(all.equal(x, 0))) {
      rlang::abort(paste0("`", nm, "` must be non-zero to determine the allometry."), call = call)
    }
    x
  }

  if (length(shape) >= 2) {
    if ("r_app_mass" %in% shape) {
      r <- nonzero(s$r_app_mass, "r_app_mass")
      b <- if ("b_sma" %in% shape) {
        s$b_sma
      } else if ("b_ols" %in% shape) {
        s$b_ols / r
      } else {
        if (isTRUE(all.equal(r, -1))) rlang::abort("`r_app_mass = -1` leaves `b_avg` uninformative.", call = call)
        2 * s$b_avg / (1 + r)
      }
    } else if (all(c("b_ols", "b_sma") %in% shape)) {
      b <- nonzero(s$b_sma, "b_sma")
      r <- s$b_ols / b
    } else if (all(c("b_ols", "b_avg") %in% shape)) {
      b <- nonzero(2 * s$b_avg - s$b_ols, "b_sma (implied by b_ols and b_avg)")
      r <- s$b_ols / b
    } else {
      b <- nonzero(s$b_sma, "b_sma")
      r <- 2 * s$b_avg / b - 1
    }
    return(c(r = r, b_sma = b))
  }

  ## One shape quantity plus both standard deviations: the SDs fix |b_sma|, the shape quantity signs it and yields r. `b_sma` alone cannot, because it is exactly what the two SDs already say.
  if (length(shape) == 1L && length(scale) == 2L) {
    ratio <- s$sd_log_append / s$sd_log_mass
    if (shape == "r_app_mass") {
      r <- nonzero(s$r_app_mass, "r_app_mass")
      return(c(r = r, b_sma = sign(r) * ratio))
    }
    if (shape == "b_ols") {
      b <- sign(nonzero(s$b_ols, "b_ols")) * ratio
      return(c(r = s$b_ols / b, b_sma = b))
    }
    if (shape == "b_avg") {
      b <- sign(nonzero(s$b_avg, "b_avg")) * ratio
      return(c(r = 2 * s$b_avg / b - 1, b_sma = b))
    }
  }

  abort_under_determined(shape, call)
}

### Reject parameter sets no joint distribution can produce, before MASS::mvrnorm fails cryptically.
check_allometry_feasible <- function(out, call) {
  if (out$sd_log_append <= 0 || out$sd_log_mass <= 0) {
    rlang::abort("Both log-scale standard deviations must be positive.", call = call)
  }
  if (abs(out$r_app_mass) > 1 + 1e-8) {
    rlang::abort(
      c(paste0("The implied `r_app_mass` is ", format(round(out$r_app_mass, 4)), ", outside [-1, 1]."),
        i = "`b_ols` and `b_sma` imply `r_app_mass = b_ols / b_sma`, so `|b_ols|` cannot exceed `|b_sma|`."),
      call = call
    )
  }
  if (isTRUE(all.equal(out$r_app_mass, 0))) {
    rlang::abort("`r_app_mass` must be non-zero: with no mass-appendage correlation the allometric slope is undefined.", call = call)
  }
  invisible(out)
}

### Anything the user pinned must agree with the solution. Recomputing and comparing catches contradictions such as b_ols/b_sma implying one correlation while r_app_mass asserts another.
check_allometry_consistent <- function(out, supplied, call) {
  bad <- character(0)
  for (nm in names(supplied)) {
    if (!isTRUE(all.equal(supplied[[nm]], out[[nm]], tolerance = 1e-6))) {
      bad <- c(bad, paste0("`", nm, "` = ", format(supplied[[nm]]),
                           ", but the other arguments imply ", format(round(out[[nm]], 6))))
    }
  }
  if (length(bad)) {
    rlang::abort(
      c("The supplied allometry is over-determined and inconsistent.", stats::setNames(bad, rep("x", length(bad)))),
      call = call
    )
  }
  invisible(out)
}

### Resolve a complete allometry from any sufficient subset of its six describing quantities.
solve_allometry <- function(r_app_mass = NULL, b_ols = NULL, b_sma = NULL, b_avg = NULL,
                            sd_log_append = NULL, sd_log_mass = NULL,
                            call = rlang::caller_env()) {
  supplied <- list(r_app_mass = r_app_mass, b_ols = b_ols, b_sma = b_sma, b_avg = b_avg,
                   sd_log_append = sd_log_append, sd_log_mass = sd_log_mass)
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]
  for (nm in names(supplied)) check_scalar_number(supplied[[nm]], nm, call)

  shape <- intersect(shape_args, names(supplied))
  scale <- intersect(scale_args, names(supplied))

  ## Shape must be pinned outright or defaulted outright: silently pairing one user-supplied slope with a default correlation would determine the other two slopes behind their back. Scale is different -- it leaves every slope and correlation untouched (it only sets raw-scale dispersion), so it takes an ordinary default.
  if (length(shape) == 0L) {
    supplied$b_sma      <- default_allometry$b_sma
    supplied$r_app_mass <- default_allometry$r_app_mass
    shape <- intersect(shape_args, names(supplied))
  }
  if (length(scale) == 0L) {
    supplied$sd_log_mass <- default_allometry$sd_log_mass
    scale <- "sd_log_mass"
  }

  rb <- solve_shape(supplied, shape, scale, call)
  r  <- unname(rb[["r"]])
  b  <- unname(rb[["b_sma"]])

  if (length(scale) == 2L) {
    sd_a <- supplied$sd_log_append
    sd_m <- supplied$sd_log_mass
  } else if ("sd_log_mass" %in% scale) {
    sd_m <- supplied$sd_log_mass
    sd_a <- abs(b) * sd_m
  } else {
    sd_a <- supplied$sd_log_append
    if (isTRUE(all.equal(b, 0))) rlang::abort("`b_sma` implied as 0; `sd_log_mass` cannot be derived.", call = call)
    sd_m <- sd_a / abs(b)
  }

  out <- list(
    r_app_mass    = r,
    b_ols         = r * b,
    b_sma         = b,
    b_avg         = b * (1 + r) / 2,
    sd_log_append = sd_a,
    sd_log_mass   = sd_m
  )
  check_allometry_feasible(out, call)
  check_allometry_consistent(out, supplied, call)
  out
}


#' Resolve an allometry from any sufficient subset of its parameters
#'
#' A log-log allometry between an appendage and body mass has three free
#' parameters: the two log-scale standard deviations and the correlation between
#' them. Six quantities describe it, and they are linked by
#'
#' \deqn{b_{SMA} = \mathrm{sign}(r)\,\sigma_A / \sigma_M, \quad
#'       b_{OLS} = r \, b_{SMA}, \quad
#'       b_{avg} = (b_{OLS} + b_{SMA}) / 2 .}
#'
#' So *any two* of `r_app_mass`, `b_ols`, `b_sma`, `b_avg` determine the other
#' two, and one standard deviation then fixes the scale.
#'
#' Supply **two shape quantities, or none**. Supplying none takes the isometric
#' default, `b_sma = 1/3, r_app_mass = 0.3`. Supplying exactly one is an error:
#' pairing your slope with a default correlation would silently determine the
#' other two slopes. Over-specification is allowed but checked for consistency.
#'
#' The standard deviations behave differently. Neither slope nor correlation
#' depends on them, so they take an ordinary default (`sd_log_mass = 0.07`) when
#' omitted. They set raw-scale dispersion and skew, which still matters, because
#' the SLI is computed on raw metrics.
#'
#' The three slopes are not interchangeable, and the difference matters:
#'
#' - `b_ols` is the slope of \eqn{E[\log A \mid \log M]}. In the bivariate
#'   lognormal that [sim_allometric()] draws from, this is the *generative*
#'   slope: the coefficient an `lm()` will recover.
#' - `b_sma` is \eqn{\sigma_A / \sigma_M}, the slope a standardised major axis
#'   fit returns. It is the exponent [calc_sli()] assumes and
#'   [build_sli_slopes_tbl()] estimates, and the one usually meant by "the
#'   allometric slope".
#' - `b_avg` is their midpoint. No estimator returns it. It is a convention for
#'   splitting the difference when OLS and SMA disagree, and is retained because
#'   earlier work used it. **Note that `b_avg = 1/3` does not mean isometry**:
#'   with `b_avg = 0.33, r_app_mass = 0.3` neither slope is near 1/3.
#'
#' @param r_app_mass Pearson correlation between `log(Append)` and `log(Mass)`.
#'   Must be non-zero.
#' @param b_ols,b_sma,b_avg The OLS, SMA, and midpoint slopes of `log(Append)`
#'   on `log(Mass)`.
#' @param sd_log_append,sd_log_mass Log-scale standard deviations.
#'
#' @return A one-row tibble with all six quantities.
#'
#' @examples
#' # Isometry as an SMA fit would report it
#' implied_allometry(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07)
#'
#' # The parameterisation used by the SLI manuscript
#' implied_allometry(b_avg = 0.33, r_app_mass = 0.3, sd_log_mass = 0.07)
#'
#' # Pin both slopes and let the correlation follow
#' implied_allometry(b_ols = 0.15, b_sma = 0.50, sd_log_mass = 0.07)
#'
#' @export
implied_allometry <- function(r_app_mass = NULL, b_ols = NULL, b_sma = NULL, b_avg = NULL,
                              sd_log_append = NULL, sd_log_mass = NULL) {
  out <- solve_allometry(
    r_app_mass = r_app_mass, b_ols = b_ols, b_sma = b_sma, b_avg = b_avg,
    sd_log_append = sd_log_append, sd_log_mass = sd_log_mass
  )
  tibble::as_tibble(out)
}
