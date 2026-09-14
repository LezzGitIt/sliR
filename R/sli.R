# Standardized Length Index ----

#' Standardized Length Index (SLI)
#'
#' Size-corrects an appendage length by rescaling every individual to a common
#' reference body mass, following the scaled mass index of Peig & Green (2009).
#' Whereas their index standardises *mass* to a reference *length*, the SLI
#' standardises a *length* to a reference *mass*:
#'
#' \deqn{SLI_i = A_i \left( \frac{M_0}{M_i} \right)^{b}}
#'
#' where \eqn{A_i} and \eqn{M_i} are individual \eqn{i}'s appendage length and
#' body mass, \eqn{M_0} is an arbitrary reference mass, and \eqn{b} is the
#' allometric scaling exponent relating length to mass.
#'
#' @param df A data frame of individual measurements.
#' @param Append,Mass Unquoted column names for appendage length and body mass.
#'   **These default to columns literally named `Append` and `Mass`.** If your
#'   columns are named otherwise, pass them explicitly (e.g. `Append = wing,
#'   Mass = mass`); an informative error is raised if the columns are absent.
#' @param b_sli Scalar allometric scaling exponent. The default `0.33`
#'   (\eqn{\approx 1/3}) assumes isometry, under which linear
#'   dimensions scale with the cube root of mass. Ignored when `control` is
#'   supplied. Set it to an empirically estimated SMA slope to relax isometry.
#' @param M0 Reference mass to which all individuals are standardised. `NULL`
#'   (the default) uses the mean of `Mass`. \eqn{M_0} is an arbitrary scalar: it
#'   rescales every SLI value by a constant factor and so changes the units of
#'   the index, never the ranking of individuals or the result of any subsequent
#'   location-invariant analysis. (Peig & Green call this \eqn{L_0}, because
#'   their index standardises to a reference *length*; the SLI inverts that.)
#' @param control Optional character vector of grouping columns (e.g.
#'   `c("Age", "Sex")`). When supplied, each individual's `b_sli` is replaced
#'   by a group-specific slope, resolved as `method` describes. Individuals
#'   whose `control` values are `NA` or match `unknown_codes` receive
#'   `sli = NA` under `method = "average"`; under `"hierarchical"` they
#'   instead cascade through the same fallback as everyone else (see
#'   [build_sli_slopes_hierarchical()]).
#' @param method How a group's slope is resolved when `control` is supplied.
#'   `"average"` (the default) uses [build_sli_slopes_tbl()]: every
#'   `control` variable's own SMA slope, averaged together unconditionally.
#'   `"hierarchical"` uses [build_sli_slopes_hierarchical()]: a
#'   reliability-gated cascade that falls back to coarser groupings (and
#'   ultimately the pooled slope) when a cell is too sparse or too weakly
#'   correlated to trust its own slope. Ignored when `control` is `NULL`.
#' @param unknown_codes Values in the `control` columns that mark an unknown
#'   group and so cannot be assigned a slope under `method = "average"`.
#' @param slope_diff_warn Passed to [build_sli_slopes_tbl()] when
#'   `method = "average"`. Warn when the per-variable SMA slopes averaged
#'   within a group differ by more than this. Only relevant with two or more
#'   `control` variables.
#' @param n_min_cell,n_min_marginal,cor_min,cor_p_max Passed to
#'   [build_sli_slopes_hierarchical()] when `method = "hierarchical"`;
#'   ignored otherwise.
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
#' @seealso [build_sli_slopes_tbl()] and [build_sli_slopes_hierarchical()] for
#'   the two ways a group's slope can be resolved, and [build_group_cor_tbl()]
#'   to check whether per-group slopes are warranted at all.
#'
#' @examples
#' set.seed(1)
#' d <- sim_allometric(n = 200)
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
#' # Per-group slopes, falling back to coarser groupings where a cell is
#' # too sparse or too weakly correlated to trust its own slope.
#' d$Sex <- rep(c("F", "M"), length.out = nrow(d))
#' head(calc_sli(d, control = "Sex", method = "hierarchical"))
#'
#' @export
calc_sli <- function(df,
                     Append = Append,
                     Mass = Mass,
                     b_sli = 0.33,
                     M0 = NULL,
                     control = NULL,
                     method = c("average", "hierarchical"),
                     unknown_codes = c("Unk", "U", "Unknown"),
                     slope_diff_warn = 0.15,
                     n_min_cell = 50,
                     n_min_marginal = 100,
                     cor_min = 0.3,
                     cor_p_max = 0.05,
                     rename_col = NULL,
                     sort = FALSE) {
  app_q   <- rlang::enquo(Append)
  mass_q  <- rlang::enquo(Mass)
  app_nm  <- rlang::as_label(app_q)
  mass_nm <- rlang::as_label(mass_q)
  check_cols(df, c(app_nm, mass_nm))
  method <- match.arg(method)

  if (!is.null(rename_col) && !rlang::is_string(rename_col)) {
    rlang::abort("`rename_col` must be a single string or `NULL`.")
  }
  if (!is.null(control) && !missing(b_sli)) {
    rlang::warn("`b_sli` is ignored when `control` is supplied; per-group SMA slopes are used instead.")
  }

  ## M0 is a scalar reference mass, not a column: injected with `!!` so a column of the same name cannot shadow it.
  if (is.null(M0)) M0 <- mean(df[[mass_nm]], na.rm = TRUE)
  if (!rlang::is_scalar_double(M0) && !rlang::is_scalar_integer(M0)) {
    rlang::abort("`M0` must be a single number or `NULL`.")
  }

  if (is.null(control)) {
    df_sli <- dplyr::mutate(df, sli = {{ Append }} * (!!M0 / {{ Mass }})^(!!b_sli))
  } else {
    slopes_tbl <- if (method == "average") {
      build_sli_slopes_tbl(
        df, Append = !!app_q, Mass = !!mass_q,
        control = control, unknown_codes = unknown_codes,
        slope_diff_warn = slope_diff_warn
      )
    } else {
      build_sli_slopes_hierarchical(
        df, Append = !!app_q, Mass = !!mass_q,
        control = control, unknown_codes = unknown_codes,
        n_min_cell = n_min_cell, n_min_marginal = n_min_marginal,
        cor_min = cor_min, cor_p_max = cor_p_max
      )
    }
    df_sli <- df |>
      dplyr::left_join(
        dplyr::select(slopes_tbl, dplyr::all_of(control), "b_sli_avg"),
        by = control
      ) |>
      dplyr::mutate(sli = {{ Append }} * (!!M0 / {{ Mass }})^.data$b_sli_avg) |>
      dplyr::select(-"b_sli_avg")
  }

  if (sort) df_sli <- dplyr::arrange(df_sli, dplyr::desc(sli))
  if (!is.null(rename_col)) df_sli <- dplyr::rename(df_sli, !!rename_col := sli)
  df_sli
}
