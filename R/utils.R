# Internal helpers ----

### Abort with a helpful message when a user-supplied (or defaulted) column name is absent from the data.
# The hint matters because `Append` and `Mass` are unquoted NSE defaults: a user whose columns are named `wing`/`mass` otherwise gets an opaque "object not found" error from deep inside a dplyr verb.
check_cols <- function(df, cols, call = rlang::caller_env()) {
  absent <- setdiff(cols, names(df))
  if (length(absent) == 0) return(invisible(df))
  header <- paste0("Column", if (length(absent) > 1) "s" else "", " not found in `df`: ",
                   paste0("`", absent, "`", collapse = ", "))
  ## The NSE hint is only relevant when the *defaulted* trait columns are what went missing; a mistyped `control` name needs no such explanation.
  hint <- if (any(c("Append", "Mass") %in% absent)) {
    c(i = "Trait columns are unquoted and default to `Append` and `Mass`.",
      i = "Supply your own names, e.g. `calc_sli(df, Append = wing, Mass = mass)`.")
  } else {
    character(0)
  }
  rlang::abort(c(header, hint), call = call)
}

### Drop rows lying more than `n_sd` standard deviations from the column mean, applied jointly across `cols`.
# Trimming happens on whatever scale `cols` are already on. On a raw (lognormal) morphometric this is asymmetric and shaves the right tail harder than the left; trim the log-scale columns instead if that matters for your use.
trim_outliers <- function(df, cols, n_sd = 3) {
  if (is.null(n_sd) || is.infinite(n_sd)) return(df)
  keep <- Reduce(`&`, lapply(cols, function(cl) {
    x <- df[[cl]]
    abs(x - mean(x, na.rm = TRUE)) <= n_sd * stats::sd(x, na.rm = TRUE)
  }))
  df[keep, , drop = FALSE]
}

### Convert a correlation matrix plus a vector of standard deviations into a covariance matrix.
# Replaces MBESS::cor2cov so the package does not take on that dependency for a one-line identity.
cor_to_cov <- function(cor_mat, sds) {
  d <- diag(sds, nrow = length(sds))
  out <- d %*% cor_mat %*% d
  dimnames(out) <- dimnames(cor_mat)
  out
}

### Reject correlation matrices that no joint distribution can realise, before MASS::mvrnorm fails cryptically.
# With a gradient the infeasible region is exactly the biologically interesting one: an appendage that lengthens while mass falls along the gradient cannot coexist with a strong positive appendage-mass correlation. Name that tension rather than reporting a bare determinant.
check_pos_def <- function(cor_mat, call = rlang::caller_env()) {
  eig <- eigen(cor_mat, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) > 0) return(invisible(cor_mat))

  hint <- if (nrow(cor_mat) == 3L) {
    grad <- colnames(cor_mat)[3]
    c(i = paste0("`r_app_mass`, `r_grad_app`, and `r_grad_mass` must satisfy ",
                 "1 - r_app_mass^2 - r_grad_app^2 - r_grad_mass^2 + ",
                 "2*r_app_mass*r_grad_app*r_grad_mass > 0."),
      i = paste0("Appendage and mass responding to ", grad,
                 " in opposite directions is incompatible with a strong positive `r_app_mass`."),
      i = "Weaken `r_app_mass`, or move `r_grad_app` and `r_grad_mass` closer together.")
  } else {
    c(i = "Shrink the correlations toward zero, or make their signs mutually consistent.")
  }

  rlang::abort(
    c("The requested correlation matrix is not positive definite.",
      i = "No joint distribution has all of these pairwise correlations simultaneously.",
      hint),
    call = call
  )
}

### Add independent measurement and transient (biological) error, each specified as a fraction of the variable's own standard deviation.
add_error <- function(x, meas_error = 0, transient_error = 0) {
  n <- length(x)
  sd_x <- stats::sd(x, na.rm = TRUE)
  x +
    stats::rnorm(n, 0, sd = sd_x * meas_error) +
    stats::rnorm(n, 0, sd = sd_x * transient_error)
}
