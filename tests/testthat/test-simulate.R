# build_cov_mat ----

test_that("build_cov_mat is symmetric, positive definite, and correctly named", {
  Sigma <- build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07)
  expect_equal(dim(Sigma), c(2L, 2L))
  expect_equal(dimnames(Sigma), list(c("Append", "Mass"), c("Append", "Mass")))
  expect_equal(Sigma, t(Sigma))
  expect_true(all(eigen(Sigma, only.values = TRUE)$values > 0))
})

test_that("build_cov_mat recovers the requested correlations and SDs", {
  Sigma <- build_cov_mat(b_sma = 0.5, r_app_mass = 0.3, sd_log_mass = 0.07)
  R <- stats::cov2cor(Sigma)
  expect_equal(R[["Append", "Mass"]], 0.3)
  expect_equal(sqrt(Sigma[["Mass", "Mass"]]), 0.07)
  expect_equal(sqrt(Sigma[["Append", "Append"]]), 0.035)
})

test_that("build_cov_mat gains a unit-SD gradient block named after the gradient", {
  Sigma <- build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07,
                         gradient = "Temperature", r_grad_app = -0.3, r_grad_mass = -0.1)
  expect_equal(dim(Sigma), c(3L, 3L))
  expect_equal(colnames(Sigma)[3], "Temperature")
  # The gradient defaults to unit SD: this matrix describes correlations, not units.
  expect_equal(Sigma[["Temperature", "Temperature"]], 1)

  R <- stats::cov2cor(Sigma)
  expect_equal(R[["Append", "Temperature"]], -0.3)
  expect_equal(R[["Mass", "Temperature"]], -0.1)
})

test_that("sd_gradient scales only the gradient block, leaving morphology and correlations", {
  base   <- build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07,
                          gradient = "Temp", r_grad_app = -0.3, r_grad_mass = -0.1)
  scaled <- build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07,
                          gradient = "Temp", r_grad_app = -0.3, r_grad_mass = -0.1,
                          sd_gradient = 0.18)
  expect_equal(sqrt(scaled[["Temp", "Temp"]]), 0.18)
  # Morphology block and every correlation are invariant to sd_gradient.
  expect_equal(base[1:2, 1:2], scaled[1:2, 1:2])
  expect_equal(stats::cov2cor(base), stats::cov2cor(scaled))
  expect_error(
    build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07,
                  gradient = "Temp", sd_gradient = -1),
    "positive number"
  )
})

test_that("build_cov_mat rejects impossible correlation triples", {
  expect_error(
    build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07,
                  gradient = "Temp", r_grad_app = 0.7, r_grad_mass = -0.7),
    "positive definite"
  )
})


# sim_allometric: allometry ----

test_that("sim_allometric returns the documented columns", {
  set.seed(3)
  d <- sim_allometric(n = 100)
  expect_named(d, c("Append", "Mass", "Append_log", "Mass_log"))
  expect_equal(d$Append_log, log(d$Append))
})

test_that("log = FALSE drops the log columns without changing the retained ones", {
  set.seed(3); with_log <- sim_allometric(n = 200, meas_error = 0.1)
  set.seed(3); no_log   <- sim_allometric(n = 200, meas_error = 0.1, log = FALSE)
  expect_named(no_log, c("Append", "Mass"))
  expect_equal(no_log$Append, with_log$Append)
  expect_equal(no_log$Mass, with_log$Mass)
})

test_that("sim_allometric hits the requested correlation and slope when untrimmed", {
  set.seed(3)
  d <- sim_allometric(n = 1000, b_sma = 1/3, r_app_mass = 0.3, trim_sd = NULL)
  expect_equal(cor(d$Append_log, d$Mass_log), 0.3, tolerance = 1e-6)
  expect_equal(unname(coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]]), 1/3, tolerance = 1e-6)
})

test_that("sim_allometric trims only when asked", {
  set.seed(3)
  expect_equal(nrow(sim_allometric(n = 500, trim_sd = NULL)), 500L)
  expect_lt(nrow(sim_allometric(n = 500, trim_sd = 2)), 500L)
})

test_that("sim_allometric is reproducible under a fixed seed", {
  set.seed(99); a <- sim_allometric(n = 200, meas_error = 0.1)
  set.seed(99); b <- sim_allometric(n = 200, meas_error = 0.1)
  expect_equal(a, b)
})

test_that("gradient arguments without a gradient name are an error, not a silent no-op", {
  expect_error(sim_allometric(n = 10, r_grad_app = -0.3), "without naming a gradient")
  expect_error(sim_allometric(n = 10, gradient_range = c(0, 1)), "without naming a gradient")
})


# sim_allometric: the gradient ----

## The four gradients a user is most likely to reach for, with the bounds each is naturally described by.
gradients <- list(
  Temperature = c(0, 24),
  Rainfall    = c(400, 2000),
  Year        = c(1970, 2026),
  Latitude    = c(-90, 90)
)

test_that("a named gradient adds exactly one correctly named column", {
  set.seed(1)
  d <- sim_allometric(n = 200, gradient = "Temperature", r_grad_app = -0.3)
  expect_named(d, c("Append", "Mass", "Temperature", "Append_log", "Mass_log"))
})

test_that("gradient must be a usable, non-colliding column name", {
  expect_error(sim_allometric(n = 10, gradient = "Mass"), "already exists")
  expect_error(sim_allometric(n = 10, gradient = ""), "non-empty string")
  expect_error(sim_allometric(n = 10, gradient = c("a", "b")), "non-empty string")
})

test_that("uniform gradients hit their target correlations exactly, at any scale", {
  for (nm in names(gradients)) {
    set.seed(3)
    d <- sim_allometric(n = 2000, gradient = nm, gradient_range = gradients[[nm]],
                        r_grad_app = -0.3, r_grad_mass = -0.1, trim_sd = NULL)
    g <- d[[nm]]
    expect_equal(cor(g, d$Append_log), -0.3, tolerance = 1e-6, info = nm)
    expect_equal(cor(g, d$Mass_log),   -0.1, tolerance = 1e-6, info = nm)
    # The gradient must not disturb the allometry it sits alongside.
    expect_equal(cor(d$Append_log, d$Mass_log), 0.3, tolerance = 1e-6, info = nm)
    expect_equal(unname(coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]]),
                 1/3, tolerance = 1e-6, info = nm)
  }
})

test_that("gradient_range fixes the mean and SD, and uniform draws stay inside the bounds", {
  for (nm in names(gradients)) {
    rng <- gradients[[nm]]
    set.seed(3)
    d <- sim_allometric(n = 2000, gradient = nm, gradient_range = rng, trim_sd = NULL)
    g <- d[[nm]]
    expect_equal(mean(g), mean(rng), tolerance = 1e-6, info = nm)
    expect_equal(sd(g), diff(rng) / sqrt(12), tolerance = 1e-6, info = nm)
    expect_gte(min(g), rng[1])
    expect_lte(max(g), rng[2])
  }
})

test_that("a uniform gradient covers its range evenly; a normal one does not", {
  set.seed(3)
  u <- sim_allometric(n = 3000, gradient = "Temperature", gradient_range = c(0, 24), trim_sd = NULL)
  tb_u <- table(cut(u$Temperature, 15))
  expect_lt(max(tb_u) / min(tb_u), 1.5)

  set.seed(3)
  n <- sim_allometric(n = 3000, gradient = "Temperature", gradient_dist = "normal",
                      mean_gradient = 12, sd_gradient = 6.93, trim_sd = NULL)
  tb_n <- table(cut(n$Temperature, 15))
  expect_gt(max(tb_n) / min(tb_n), 100)
})

test_that("a normal gradient is unbounded, which is why uniform is the default", {
  set.seed(3)
  d <- sim_allometric(n = 3000, gradient = "Latitude", gradient_dist = "normal",
                      mean_gradient = 0, sd_gradient = 51.96, trim_sd = NULL)
  expect_gt(mean(abs(d$Latitude) > 90), 0.01)
})

test_that("gradient_dist = 'normal' reproduces the target covariance", {
  set.seed(3)
  d <- sim_allometric(n = 3000, gradient = "Temperature", gradient_dist = "normal",
                      r_grad_app = -0.3, r_grad_mass = -0.1, trim_sd = NULL)
  expect_equal(cor(d$Temperature, d$Append_log), -0.3, tolerance = 1e-6)
  expect_equal(cor(d$Temperature, d$Mass_log),   -0.1, tolerance = 1e-6)
  expect_equal(cor(d$Append_log, d$Mass_log),     0.3, tolerance = 1e-6)

  target <- stats::cov2cor(build_cov_mat(gradient = "Temperature",
                                         r_grad_app = -0.3, r_grad_mass = -0.1))
  got <- stats::cor(cbind(Append = d$Append_log, Mass = d$Mass_log, Temperature = d$Temperature))
  expect_equal(got, target, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("empirical = FALSE relaxes exactness but keeps the structure", {
  set.seed(3)
  d <- sim_allometric(n = 4000, gradient = "Temperature", gradient_range = c(0, 24),
                      r_grad_app = -0.3, trim_sd = NULL, empirical = FALSE)
  expect_equal(cor(d$Temperature, d$Append_log), -0.3, tolerance = 0.05)
  expect_gte(min(d$Temperature), 0)
  expect_lte(max(d$Temperature), 24)
})

test_that("gradient_range is refused where it makes no sense", {
  expect_error(sim_allometric(n = 10, gradient = "Year", gradient_dist = "normal",
                              gradient_range = c(1970, 2026)), "requires .*uniform")
  expect_error(sim_allometric(n = 10, gradient = "Year", gradient_range = c(1970, 2026),
                              sd_gradient = 5), "not both")
  expect_error(sim_allometric(n = 10, gradient = "Year", gradient_range = c(2026, 1970)),
               "two increasing numbers")
})

test_that("impossible gradient correlations are rejected with a biological explanation", {
  expect_error(
    sim_allometric(n = 100, gradient = "Temperature", r_app_mass = 0.3, b_sma = 1/3,
                   sd_log_mass = 0.07, r_grad_app = 0.7, r_grad_mass = -0.7),
    "opposite directions"
  )
})

test_that("trimming drops rows from the gradient column too, keeping the tibble rectangular", {
  set.seed(3)
  d <- sim_allometric(n = 500, gradient = "Temperature", trim_sd = 2)
  expect_lt(nrow(d), 500L)
  expect_false(anyNA(d$Temperature))
})


# sim_grid ----

test_that("sim_grid runs one scenario per row of params and stamps .scenario", {
  grid <- expand.grid(n = 50, r_app_mass = c(0.1, 0.3, 0.5), b_sma = 1/3)
  set.seed(1)
  sims <- sim_grid(grid)
  expect_named(sims, c(".scenario", "Append", "Mass", "Append_log", "Mass_log"))
  expect_equal(unique(sims$.scenario), 1:3)
  expect_true(is.integer(sims$.scenario))
  expect_equal(as.integer(table(sims$.scenario)), rep(50L, 3))
})

test_that("sim_grid forwards constant arguments via ... to every scenario", {
  grad_grid <- expand.grid(n = 50, r_grad_app = c(-0.3, 0, 0.3), r_grad_mass = -0.1)
  set.seed(1)
  sims <- sim_grid(grad_grid, gradient = "Temperature", gradient_range = c(0, 24), trim_sd = NULL)
  expect_true("Temperature" %in% names(sims))
  expect_equal(as.integer(table(sims$.scenario)), rep(50L, 3))
  expect_gte(min(sims$Temperature), 0)
  expect_lte(max(sims$Temperature), 24)
})

test_that("sim_grid matches calling sim_allometric row by row", {
  grid <- expand.grid(n = 30, r_app_mass = c(0.2, 0.4), b_sma = 1/3)
  set.seed(7); combined <- sim_grid(grid)
  set.seed(7)
  separate <- dplyr::bind_rows(purrr::pmap(grid, sim_allometric), .id = "x")
  expect_equal(dplyr::select(combined, -".scenario"), dplyr::select(separate, -"x"))
})

test_that("sim_grid rejects a bad params or id_col", {
  expect_error(sim_grid(data.frame()), "at least one row")
  expect_error(sim_grid(list(n = 50)), "data frame")
  grid <- expand.grid(n = 10, r_app_mass = 0.3, b_sma = 1/3)
  expect_error(sim_grid(grid, id_col = "n"), "already names a column")
  expect_error(sim_grid(grid, id_col = c("a", "b")), "single non-empty string")
})


# sim_correlated ----

test_that("sim_correlated hits the requested correlation and moments with empirical = TRUE", {
  set.seed(5)
  d <- sim_correlated(n = 1000, r_app_mass = 0.3, mu_append = 180, mu_mass = 80,
                      sd_append = 10, sd_mass = 5)
  expect_named(d, c("Append", "Mass"))
  expect_equal(cor(d$Append, d$Mass), 0.3, tolerance = 1e-8)
  expect_equal(mean(d$Append), 180, tolerance = 1e-8)
  expect_equal(sd(d$Mass), 5, tolerance = 1e-8)
})

test_that("transient error in mass alone attenuates the mass-appendage correlation", {
  set.seed(5); clean <- sim_correlated(n = 2000, r_app_mass = 0.3)
  set.seed(5); noisy <- sim_correlated(n = 2000, r_app_mass = 0.3, transient_error_mass = 1)
  expect_lt(abs(cor(noisy$Append, noisy$Mass)), abs(cor(clean$Append, clean$Mass)))
})

test_that("the deprecated `r` argument still works and warns", {
  ## `.frequency = "once"` fires at most once per session, so this must be the only
  ## call to `sim_correlated(r = ...)` in the test suite; a second call would not warn.
  set.seed(5)
  expect_warning(via_r <- sim_correlated(n = 500, r = 0.3), "deprecated")
  set.seed(5)
  via_r_app_mass <- sim_correlated(n = 500, r_app_mass = 0.3)
  expect_equal(via_r, via_r_app_mass)
})
