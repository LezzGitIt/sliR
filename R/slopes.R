# Per-group allometric slopes ----

### Slope tables are keyed on the control columns and later joined back onto the individual-level data, so the key's type must survive the round trip through smatr, which returns group levels as row names (always character).
restore_key_type <- function(x, template) {
  if (is.factor(template)) factor(x, levels = levels(template)) else x
}

### smatr fits a separate slope per level only for factor/character groups; a numeric control column would be treated as a covariate and silently return a single slope.
check_control_types <- function(df, control, call = rlang::caller_env()) {
  bad <- control[!vapply(df[control], \(x) is.character(x) || is.factor(x), logical(1))]
  if (length(bad)) {
    rlang::abort(
      c(paste0("`control` columns must be character or factor: ",
               paste0("`", bad, "`", collapse = ", "), " ",
               if (length(bad) > 1) "are" else "is", " not."),
        i = "A numeric grouping column would be fitted as a covariate, yielding one slope rather than one per group."),
      call = call
    )
  }
  invisible(df)
}

#' Per-group SMA allometric slopes
#'
#' Fits one standardised major axis (SMA) regression of `log(Append)` on
#' `log(Mass)` per `control` variable, letting the slope vary across that
#' variable's levels. Every combination of levels present in the data then
#' receives `b_sli_avg`, the mean of the slopes its levels attract.
#'
#' Note the structure: with `control = c("Age", "Sex")` this fits **two**
#' models — one interacting mass with age, one interacting mass with sex — and
#' averages their slopes per age × sex cell. It does not fit a single
#' age × sex interaction. That keeps the per-cell sample sizes at the marginal
#' rather than the joint level, which matters when some cells are sparse.
#'
#' Rows whose `control` values are `NA` or match `unknown_codes` are excluded
#' from slope fitting and from the returned table.
#'
#' @inheritParams calc_sli
#' @param control Character vector of grouping columns. Each must be character
#'   or factor.
#'
#' @return A tibble with one row per observed combination of `control` levels:
#'   the control columns, `n` (rows contributing to that cell), one
#'   `b_sma_<var>` column per control variable, and their row mean
#'   `b_sli_avg`. A cell whose slope is missing for any control variable gets
#'   `b_sli_avg = NA`.
#'
#' @seealso [build_group_cor_tbl()], which flags groups whose mass–appendage
#'   relationship is too weak for a per-group slope to be meaningful.
#'
#' @examples
#' set.seed(1)
#' d <- sim_allometric(n = 400)
#' d$Sex <- rep(c("F", "M"), length.out = nrow(d))
#' build_sli_slopes_tbl(d, control = "Sex")
#'
#' @export
build_sli_slopes_tbl <- function(df,
                                 Append = Append,
                                 Mass = Mass,
                                 control,
                                 unknown_codes = c("Unk", "U", "Unknown")) {
  app_nm  <- rlang::as_label(rlang::enquo(Append))
  mass_nm <- rlang::as_label(rlang::enquo(Mass))
  check_cols(df, c(app_nm, mass_nm, control))
  check_control_types(df, control)

  df_fit <- df |>
    dplyr::mutate(
      .log_app  = log(.data[[app_nm]]),
      .log_mass = log(.data[[mass_nm]])
    ) |>
    dplyr::filter(dplyr::if_all(
      dplyr::all_of(control),
      \(x) !is.na(x) & !as.character(x) %in% unknown_codes
    ))

  ## One SMA per control variable, collecting that variable's per-level slopes.
  slope_tbls <- purrr::map(control, \(grp_var) {
    fmla <- stats::as.formula(paste(".log_app ~ .log_mass *", grp_var))
    fit  <- smatr::sma(fmla, data = df_fit, method = "SMA")
    tbl  <- tibble::as_tibble(stats::coef(fit), rownames = grp_var)
    tbl[[grp_var]] <- restore_key_type(tbl[[grp_var]], df_fit[[grp_var]])
    tbl |>
      dplyr::rename(!!paste0("b_sma_", grp_var) := "slope") |>
      dplyr::select(dplyr::all_of(grp_var), dplyr::starts_with("b_sma_"))
  })

  ## Every level combination observed in the known data, with its sample size.
  slopes_tbl <- dplyr::count(df_fit, !!!rlang::syms(control))
  for (i in seq_along(control)) {
    slopes_tbl <- dplyr::left_join(slopes_tbl, slope_tbls[[i]], by = control[[i]])
  }

  slope_cols <- paste0("b_sma_", control)
  dplyr::mutate(
    slopes_tbl,
    b_sli_avg = rowMeans(dplyr::across(dplyr::all_of(slope_cols)), na.rm = FALSE)
  )
}


#' Per-group mass–appendage correlation diagnostic
#'
#' Reports, for each combination of `control` levels, the OLS regression of
#' mass on appendage length together with its correlation and p-value. Use it
#' before [calc_sli()]`(control = ...)`: a group whose mass and appendage
#' length are uncorrelated has no meaningful allometric relationship, so the
#' SMA slope fitted to it is noise and should not be used to size-correct its
#' members.
#'
#' @inheritParams build_sli_slopes_tbl
#'
#' @return A tibble with one row per observed combination of `control` levels:
#'   the control columns, `n`, `b_ols` (slope of mass on appendage), `r`
#'   (Pearson correlation), and `p_value`.
#'
#' @examples
#' set.seed(1)
#' d <- sim_allometric(n = 400)
#' d$Sex <- rep(c("F", "M"), length.out = nrow(d))
#' build_group_cor_tbl(d, control = "Sex")
#'
#' @export
build_group_cor_tbl <- function(df,
                                Append = Append,
                                Mass = Mass,
                                control,
                                unknown_codes = c("Unk", "U", "Unknown")) {
  app_nm  <- rlang::as_label(rlang::enquo(Append))
  mass_nm <- rlang::as_label(rlang::enquo(Mass))
  check_cols(df, c(app_nm, mass_nm, control))

  df |>
    dplyr::filter(dplyr::if_all(
      dplyr::all_of(control),
      \(x) !is.na(x) & !as.character(x) %in% unknown_codes
    )) |>
    dplyr::group_by(!!!rlang::syms(control)) |>
    dplyr::group_modify(\(grp, key) {
      fmla <- stats::as.formula(paste("`", mass_nm, "` ~ `", app_nm, "`", sep = ""))
      mod  <- stats::lm(fmla, data = grp)
      tibble::tibble(
        n       = nrow(grp),
        b_ols   = stats::coef(mod)[[2]],
        r       = stats::cor(grp[[app_nm]], grp[[mass_nm]], use = "complete.obs"),
        p_value = summary(mod)$coefficients[2, "Pr(>|t|)"]
      )
    }) |>
    dplyr::ungroup()
}
