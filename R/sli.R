# Standardized Length Index ----

#' Standardized Length Index (SLI)
#'
#' Size-corrects an appendage length by rescaling every individual to a common
#' reference body mass, following the scaled mass index of Peig & Green (2009).
#' Whereas their index standardises *mass* to a reference *length*, the SLI
#' standardises a *length* to a reference *mass*:
#'
#' \deqn{SLI_i = A_i \left( \frac{L_0}{M_i} \right)^{b}}
#'
#' where \eqn{A_i} and \eqn{M_i} are individual \eqn{i}'s appendage length and
#' body mass, \eqn{L_0} is an arbitrary reference mass, and \eqn{b} is the
#' allometric scaling exponent relating length to mass.
#'
#' @param df A data frame of individual measurements.
#' @param Append,Mass Unquoted column names for appendage length and body mass.
#'   **These default to columns literally named `Append` and `Mass`.** If your
#'   columns are named otherwise, pass them explicitly (e.g. `Append = wing,
#'   Mass = mass`); an informative error is raised if the columns are absent.
#' @param b_sli Scalar allometric scaling exponent. The default `0.33`
#'   (\eqn{\approx 1/3}) assumes geometric similarity, under which linear
#'   dimensions scale with the cube root of mass. Ignored when `control` is
#'   supplied. Set it to an empirically estimated SMA slope to relax isometry.
#' @param L0 Reference mass to which all individuals are standardised. `NULL`
#'   (the default) uses the mean of `Mass`. Peig & Green treat \eqn{L_0} as an
#'   arbitrary scalar: it rescales every SLI value by a constant factor and so
#'   changes the units of the index, never the ranking of individuals or the
#'   result of any subsequent location-invariant analysis.
#' @param control Optional character vector of grouping columns (e.g.
#'   `c("Age", "Sex")`). When supplied, `build_sli_slopes_tbl()` estimates a
#'   separate SMA slope per group and each individual receives their group's
#'   averaged slope in place of `b_sli`. Individuals whose `control` values are
#'   `NA` or match `unknown_codes` receive `sli = NA`.
#' @param unknown_codes Values in the `control` columns that mark an unknown
#'   group and so cannot be assigned a slope.
#' @param rename_col Optional string naming the output column. `NULL` (the
#'   default) leaves it as `sli`.
#' @param sort If `TRUE`, return rows sorted by descending SLI. Defaults to
#'   `FALSE`, preserving the input row order.
#'
#' @return `df` with one additional numeric column (`sli`, or `rename_col`).
#'
#' @references
#' Peig, J. & Green, A.J. (2009) New perspectives for estimating body condition
#' from mass/length data: the scaled mass index as an alternative method.
#' *Oikos* 118, 1883–1891. \doi{10.1111/j.1600-0706.2009.17643.x}
#'
#' @seealso [build_sli_slopes_tbl()] for the per-group slopes, and
#'   [build_group_cor_tbl()] to check whether per-group slopes are warranted.
#'
#' @examples
#' set.seed(1)
#' d <- sim_allometric(n = 200, b_avg = 0.33)
#'
#' # Isometric standardisation
#' head(calc_sli(d, b_sli = 0.33))
#'
#' # Two indices side by side, in the original row order
#' d |>
#'   calc_sli(b_sli = 0.33, rename_col = "sli_isometry") |>
#'   calc_sli(b_sli = 0.25, rename_col = "sli_shallow") |>
#'   head()
#'
#' @export
calc_sli <- function(df,
                     Append = Append,
                     Mass = Mass,
                     b_sli = 0.33,
                     L0 = NULL,
                     control = NULL,
                     unknown_codes = c("Unk", "U", "Unknown"),
                     rename_col = NULL,
                     sort = FALSE) {
  app_q   <- rlang::enquo(Append)
  mass_q  <- rlang::enquo(Mass)
  app_nm  <- rlang::as_label(app_q)
  mass_nm <- rlang::as_label(mass_q)
  check_cols(df, c(app_nm, mass_nm))

  if (!is.null(rename_col) && !rlang::is_string(rename_col)) {
    rlang::abort("`rename_col` must be a single string or `NULL`.")
  }
  if (!is.null(control) && !missing(b_sli)) {
    rlang::warn("`b_sli` is ignored when `control` is supplied; per-group SMA slopes are used instead.")
  }

  ## L0 is a scalar reference mass, not a column: injected with `!!` so a column of the same name cannot shadow it.
  if (is.null(L0)) L0 <- mean(df[[mass_nm]], na.rm = TRUE)
  if (!rlang::is_scalar_double(L0) && !rlang::is_scalar_integer(L0)) {
    rlang::abort("`L0` must be a single number or `NULL`.")
  }

  if (is.null(control)) {
    df_sli <- dplyr::mutate(df, sli = {{ Append }} * (!!L0 / {{ Mass }})^(!!b_sli))
  } else {
    slopes_tbl <- build_sli_slopes_tbl(
      df, Append = !!app_q, Mass = !!mass_q,
      control = control, unknown_codes = unknown_codes
    )
    df_sli <- df |>
      dplyr::left_join(
        dplyr::select(slopes_tbl, dplyr::all_of(control), "b_sli_avg"),
        by = control
      ) |>
      dplyr::mutate(sli = {{ Append }} * (!!L0 / {{ Mass }})^.data$b_sli_avg) |>
      dplyr::select(-"b_sli_avg")
  }

  if (sort) df_sli <- dplyr::arrange(df_sli, dplyr::desc(sli))
  if (!is.null(rename_col)) df_sli <- dplyr::rename(df_sli, !!rename_col := sli)
  df_sli
}


#' Ratio of squared coefficients of variation
#'
#' Computes \eqn{\lambda = CV(y)^2 / CV(x)^2}, the empirical coefficient
#' relating the variance of two traits on the log scale. With `x` a body mass
#' and `y` an appendage length, \eqn{\sqrt{\lambda}} is the standardised major
#' axis (SMA) slope, so `calc_lambda()` gives a quick moment-based check on a
#' fitted allometric exponent without fitting a model.
#'
#' @param x,y Numeric vectors of equal length. Both must be strictly positive
#'   for the coefficient of variation to be meaningful.
#' @param na.rm Drop missing values before computing means and standard
#'   deviations.
#'
#' @return A single numeric value.
#'
#' @examples
#' set.seed(1)
#' d <- sim_allometric(n = 500, b_avg = 0.33)
#' calc_lambda(x = d$Mass, y = d$Append)
#'
#' @export
calc_lambda <- function(x, y, na.rm = FALSE) {
  cv_x <- stats::sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
  cv_y <- stats::sd(y, na.rm = na.rm) / mean(y, na.rm = na.rm)
  (cv_y^2) / (cv_x^2)
}
