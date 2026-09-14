# Per-group allometric slopes ----

### Slope tables are keyed on the control columns and later joined back onto the individual-level data, so the key's type must survive the round trip through smatr, which returns group levels as row names (always character).
restore_key_type <- function(x, template) {
  if (is.factor(template)) factor(x, levels = levels(template)) else x
}

### Collapses NA and unknown_codes to one explicit "Unk" class, rather than dropping those rows. Used by build_sli_slopes_hierarchical() so unknown-coded individuals still cascade through the same levels as everyone else, instead of being silently excluded (build_sli_slopes_tbl()/build_group_cor_tbl() drop them instead, via the same unknown_codes argument used differently).
collapse_unknown_class <- function(x, unknown_codes) {
  x_chr <- as.character(x)
  dplyr::if_else(is.na(x_chr) | x_chr %in% unknown_codes, "Unk", x_chr)
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
#' age × sex interaction (which is not possible in the smatr package). That
#' keeps the per-cell sample sizes at the marginal rather than the joint level,
#' which matters when some cells are sparse.
#'
#' The cost of that averaging is that a cell whose control variables disagree
#' about the slope is collapsed to their midpoint, and nothing in `b_sli_avg`
#' records the disagreement. `slope_diff_warn` guards against this: when any
#' cell's per-variable slopes span more than that much, a warning names the
#' worst cell and its spread. Inspect the `b_sma_<var>` columns before trusting
#' `b_sli_avg`, and consider correcting on one control variable at a time.
#'
#' Rows whose `control` values are `NA` or match `unknown_codes` are excluded
#' from slope fitting and from the returned table.
#'
#' @inheritParams calc_sli
#' @param control Character vector of grouping columns. Each must be character
#'   or factor.
#' @param slope_diff_warn Warn when the per-variable SMA slopes being averaged
#'   within a cell span more than this. Defaults to `0.15`, roughly half the
#'   isometric exponent. `NULL` or `Inf` disables the check. Has no effect with
#'   a single `control` variable, where nothing is averaged.
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
                                 unknown_codes = c("Unk", "U", "Unknown"),
                                 slope_diff_warn = 0.15) {
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
  slopes_tbl <- dplyr::mutate(
    slopes_tbl,
    b_sli_avg = rowMeans(dplyr::across(dplyr::all_of(slope_cols)), na.rm = FALSE)
  )
  warn_slope_spread(slopes_tbl, control, slope_cols, slope_diff_warn)
  slopes_tbl
}


### Averaging hides disagreement: a cell whose Age slope is 0.30 and Sex slope 0.60 gets b_sli_avg = 0.45, a value neither control variable supports. Surface the widest such cell rather than letting it pass silently.
warn_slope_spread <- function(slopes_tbl, control, slope_cols, tol) {
  if (length(control) < 2 || is.null(tol) || is.infinite(tol)) return(invisible(NULL))

  slope_mat <- as.matrix(slopes_tbl[slope_cols])
  spread    <- apply(slope_mat, 1, \(x) if (anyNA(x)) NA_real_ else max(x) - min(x))
  if (all(is.na(spread)) || max(spread, na.rm = TRUE) <= tol) return(invisible(NULL))

  n_bad <- sum(spread > tol, na.rm = TRUE)
  worst <- which.max(spread)
  cell  <- paste(
    purrr::map_chr(control, \(v) paste0(v, " = ", as.character(slopes_tbl[[v]][worst]))),
    collapse = ", "
  )
  slopes <- paste(
    purrr::map_chr(slope_cols, \(cl) paste0(cl, " = ", format(round(slopes_tbl[[cl]][worst], 3), nsmall = 3))),
    collapse = ", "
  )

  rlang::warn(c(
    paste0(n_bad, " of ", nrow(slopes_tbl), " group",
           if (nrow(slopes_tbl) > 1) "s" else "",
           " average SMA slopes that differ by more than ", tol, "."),
    i = paste0("Widest (", cell, "): ", slopes,
               "; spread = ", format(round(spread[worst], 3), nsmall = 3),
               ", averaged to b_sli_avg = ",
               format(round(slopes_tbl$b_sli_avg[worst], 3), nsmall = 3), "."),
    i = "`b_sli_avg` is a midpoint neither control variable supports. Inspect the `b_sma_*` columns, or correct on one control variable at a time.",
    i = "Raise `slope_diff_warn`, or set it to `NULL`, to silence this."
  ))
  invisible(NULL)
}


#' Per-group mass–appendage correlation diagnostic
#'
#' Reports, for each combination of `control` levels, the OLS regression of
#' mass on appendage length together with its correlation and p-value. Use it
#' before [calc_sli()]`(control = ...)`: a group whose mass and appendage
#' length are uncorrelated has no meaningful SMA slope, and should not be used to
#' size-correct its members.
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


### One row of n / r / p_value / slope for a single subset -- the unit both cascade levels below are built from. r/p_value on the raw scale (matching build_group_cor_tbl()'s convention); slope is the log-scale SMA fit calc_sli() actually uses. Returns NA slope/r/p_value, not an error, when a subset is too small to fit (n < 3) -- callers decide via `reliable()` whether NA propagates further, not this helper.
fit_slope_reliability <- function(sub, app_nm, mass_nm) {
  if (nrow(sub) < 3) {
    return(tibble::tibble(n = nrow(sub), r = NA_real_, p_value = NA_real_, slope = NA_real_))
  }
  r <- suppressWarnings(stats::cor(sub[[app_nm]], sub[[mass_nm]]))
  p_value <- tryCatch(
    stats::cor.test(sub[[app_nm]], sub[[mass_nm]])$p.value,
    error = \(e) NA_real_
  )
  fmla  <- stats::as.formula(paste0(".log_app ~ .log_mass"))
  slope <- tryCatch(
    unname(stats::coef(smatr::sma(fmla, data = sub, method = "SMA"))["slope"]),
    error = \(e) NA_real_
  )
  tibble::tibble(n = nrow(sub), r = r, p_value = p_value, slope = slope)
}

### TRUE where a subset clears both a minimum sample size and a reliable mass~appendage correlation. Vectorised over a table's rows, not a single subset.
reliable <- function(tbl, n_min, cor_min, cor_p_max) {
  !is.na(tbl$r) & tbl$n >= n_min & tbl$r >= cor_min & !is.na(tbl$p_value) & tbl$p_value < cor_p_max
}


#' Per-cell SMA allometric slopes, with reliability-gated fallback
#'
#' An alternative to [build_sli_slopes_tbl()] for when some `control` cells
#' are too sparse or too weakly correlated to fit their own SMA slope. Rather
#' than averaging every `control` variable's marginal slope unconditionally,
#' this resolves one slope per cell from a cascade, finest to coarsest,
#' stopping at the first level that clears both a minimum sample size and a
#' reliable mass-appendage correlation (`cor_min`, `cor_p_max`):
#'
#' 1. **The cell itself** -- the full cross of both `control` variables (only
#'    possible, and only attempted, with exactly two). Gated by `n_min_cell`.
#' 2. **Single-variable marginals** -- each `control` variable's own slope,
#'    pooling across the other. Gated by `n_min_marginal`; averaged together
#'    when both pass, used unblended when only one does.
#' 3. **The pooled slope** -- fit on every row of `df`, regardless of `control`.
#'    Reached when no finer level passes, and used for every cell when
#'    `control` is empty.
#'
#' `df`'s pooled mass-appendage correlation gates the whole cascade: when it
#' is not itself reliable, every cell's slope is `NA` (`pass_level = "none"`)
#' rather than falling back to an unreliable pooled fit. This keeps a caller's
#' downstream SLI defined for either all of `df`'s rows or none, matching the
#' shared-N convention other body-size-correction methods need for a fair
#' comparison across them.
#'
#' Unlike [build_sli_slopes_tbl()] and [build_group_cor_tbl()], rows whose
#' `control` values are `NA` or match `unknown_codes` are **not** dropped:
#' they are collapsed to an explicit `"Unk"` class per variable, which
#' cascades through the same three levels as every other class -- not a
#' special case or an automatic fallback to pooled.
#'
#' @inheritParams build_sli_slopes_tbl
#' @param control Character vector of **0, 1, or 2** grouping columns. Each
#'   must be character or factor. Level 1 (the combined cell) is only
#'   attempted with exactly two; with `character(0)`, every row resolves to
#'   the pooled slope.
#' @param n_min_cell Minimum `n` for level 1 (the combined cell). Defaults to
#'   `50`.
#' @param n_min_marginal Minimum `n` for level 2 (a single-variable marginal).
#'   Defaults to `100`.
#' @param cor_min,cor_p_max A level passes only when its raw-scale Pearson
#'   correlation is at least `cor_min` (default `0.3`) with `p < cor_p_max`
#'   (default `0.05`).
#'
#' @return A tibble with one row per observed combination of `control`'s
#'   collapsed classes (or one row overall, when `control` is empty): the
#'   control columns, `pass_level` (`"cell"`, `"marginal"`, `"pooled"`, or
#'   `"none"`), and `b_sli_avg` -- the resolved slope, named to match
#'   [build_sli_slopes_tbl()]'s output so both feed [calc_sli()]`(control =
#'   ...)` identically. `NA` wherever `pass_level` is `"none"`.
#'
#' @seealso [build_sli_slopes_tbl()], the unconditional-averaging alternative
#'   this generalises; [build_group_cor_tbl()] for the same per-cell
#'   reliability diagnostic without the cascade.
#'
#' @examples
#' set.seed(1)
#' d <- sim_allometric(n = 500)
#' d$Age <- sample(c("Juvenile", "Adult", NA), nrow(d), replace = TRUE, prob = c(.1, .7, .2))
#' d$Sex <- sample(c("F", "M"), nrow(d), replace = TRUE)
#' build_sli_slopes_hierarchical(d, control = c("Age", "Sex"))
#'
#' @export
build_sli_slopes_hierarchical <- function(df,
                                          Append = Append,
                                          Mass = Mass,
                                          control = character(0),
                                          n_min_cell = 50,
                                          n_min_marginal = 100,
                                          cor_min = 0.3,
                                          cor_p_max = 0.05,
                                          unknown_codes = c("Unk", "U", "Unknown")) {
  app_nm  <- rlang::as_label(rlang::enquo(Append))
  mass_nm <- rlang::as_label(rlang::enquo(Mass))
  check_cols(df, c(app_nm, mass_nm, control))
  if (length(control) > 0) check_control_types(df, control)
  if (length(control) > 2) {
    rlang::abort('`control` must have length 0, 1, or 2: the combined cell (level 1) is only defined for two variables.')
  }

  d <- df |>
    dplyr::mutate(.log_app = base::log(.data[[app_nm]]), .log_mass = base::log(.data[[mass_nm]]))

  ## Level 3 / species-level gate: the pooled fit, and the final fallback for every cell.
  pooled <- fit_slope_reliability(d, app_nm, mass_nm)
  if (!isTRUE(reliable(pooled, n_min = 0, cor_min, cor_p_max))) {
    if (length(control) == 0) {
      return(tibble::tibble(pass_level = "none", b_sli_avg = NA_real_))
    }
    cells <- dplyr::count(d, !!!rlang::syms(control))[control]
    return(dplyr::mutate(cells, pass_level = "none", b_sli_avg = NA_real_))
  }
  if (length(control) == 0) {
    return(tibble::tibble(pass_level = "pooled", b_sli_avg = pooled$slope))
  }

  for (v in control) d[[paste0(".", v, "_cls")]] <- collapse_unknown_class(d[[v]], unknown_codes)
  cls_cols <- paste0(".", control, "_cls")

  ## Level 1: the combined cell, only defined with two control variables.
  if (length(control) == 2) {
    cell_tbl <- d |>
      dplyr::group_by(dplyr::across(dplyr::all_of(cls_cols))) |>
      dplyr::group_modify(\(grp, key) fit_slope_reliability(grp, app_nm, mass_nm)) |>
      dplyr::ungroup()
    cell_tbl$l1_pass  <- reliable(cell_tbl, n_min_cell, cor_min, cor_p_max)
    cell_tbl$l1_slope <- cell_tbl$slope
    d <- dplyr::left_join(d, dplyr::select(cell_tbl, dplyr::all_of(cls_cols), "l1_slope", "l1_pass"), by = cls_cols)
  } else {
    d$l1_slope <- NA_real_
    d$l1_pass  <- FALSE
  }

  ## Level 2: one marginal fit per control variable, masked to NA wherever it fails -- so
  ## a per-row average across the masked columns automatically reduces to whichever
  ## marginal(s) passed, unblended when only one does.
  for (v in control) {
    cls <- paste0(".", v, "_cls")
    marg_tbl <- d |>
      dplyr::group_by(dplyr::across(dplyr::all_of(cls))) |>
      dplyr::group_modify(\(grp, key) fit_slope_reliability(grp, app_nm, mass_nm)) |>
      dplyr::ungroup()
    marg_tbl$pass <- reliable(marg_tbl, n_min_marginal, cor_min, cor_p_max)
    marg_tbl[[paste0(".", v, "_marg_masked")]] <- dplyr::if_else(marg_tbl$pass, marg_tbl$slope, NA_real_)
    d <- dplyr::left_join(d, dplyr::select(marg_tbl, dplyr::all_of(cls), dplyr::ends_with("_marg_masked")), by = cls)
  }

  masked_cols <- paste0(".", control, "_marg_masked")
  marg_avg    <- rowMeans(as.data.frame(d[masked_cols]), na.rm = TRUE)  # NaN where none passed
  n_marg_pass <- rowSums(!is.na(as.data.frame(d[masked_cols])))

  d$b_sli_avg   <- dplyr::case_when(
    d$l1_pass %in% TRUE ~ d$l1_slope,
    n_marg_pass >= 1     ~ marg_avg,
    TRUE                 ~ pooled$slope
  )
  d$pass_level <- dplyr::case_when(
    d$l1_pass %in% TRUE ~ "cell",
    n_marg_pass >= 1     ~ "marginal",
    TRUE                 ~ "pooled"
  )

  d |>
    dplyr::distinct(dplyr::across(dplyr::all_of(control)), .data$pass_level, .data$b_sli_avg) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(control)))
}
