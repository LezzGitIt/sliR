test_that("build_cov_mat is symmetric, positive definite, and correctly named", {
  Sigma <- build_cov_mat(b_avg = 0.33, r_am = 0.5, r_at = -0.1, r_mt = -0.4)
  expect_equal(dim(Sigma), c(3L, 3L))
  expect_equal(dimnames(Sigma), list(c("Append", "Mass", "Temp"), c("Append", "Mass", "Temp")))
  expect_equal(Sigma, t(Sigma))
  expect_true(all(eigen(Sigma, only.values = TRUE)$values > 0))
})

test_that("build_cov_mat recovers the requested correlations", {
  Sigma <- build_cov_mat(b_avg = 0.33, r_am = 0.5, r_at = -0.1, r_mt = -0.4)
  R     <- stats::cov2cor(Sigma)
  expect_equal(R[["Append", "Mass"]], 0.5)
  expect_equal(R[["Append", "Temp"]], -0.1)
  expect_equal(R[["Mass", "Temp"]], -0.4)
})

test_that("build_cov_mat places sd_log_morph on the trait named by vary", {
  by_mass <- build_cov_mat(sd_log_morph = 0.07, vary = "mass")
  expect_equal(sqrt(by_mass[["Mass", "Mass"]]), 0.07)

  by_app <- build_cov_mat(sd_log_morph = 0.07, vary = "append")
  expect_equal(sqrt(by_app[["Append", "Append"]]), 0.07)
})

test_that("build_cov_mat rejects degenerate and impossible parameters", {
  expect_error(build_cov_mat(r_am = 0), "non-zero")
  expect_error(build_cov_mat(b_avg = 0, vary = "append"), "non-zero")
  expect_error(build_cov_mat(r_am = 0.9, r_at = 0.9, r_mt = -0.9), "positive definite")
})

test_that("sim_allometric hits the requested b_avg, the mean of the OLS and SMA slopes", {
  set.seed(3)
  d <- sim_allometric(n = 2000, b_avg = 0.33, r_am = 0.3, trim_sd = NULL)
  b_ols <- unname(coef(stats::lm(Append_log ~ Mass_log, data = d))[2])
  b_sma <- unname(coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]])
  expect_equal(mean(c(b_ols, b_sma)), 0.33, tolerance = 1e-6)
})

test_that("sim_allometric hits the requested correlation when untrimmed", {
  set.seed(3)
  d <- sim_allometric(n = 1000, r_am = 0.3, trim_sd = NULL)
  expect_equal(cor(d$Append_log, d$Mass_log), 0.3, tolerance = 1e-6)
})

test_that("sim_allometric returns the documented columns, logged after error is added", {
  set.seed(3)
  d <- sim_allometric(n = 100, meas_error = 0.2)
  expect_named(d, c("Append", "Mass", "Temp", "Append_log", "Mass_log"))
  expect_equal(d$Append_log, log(d$Append))
  expect_equal(d$Mass_log, log(d$Mass))
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

test_that("sim_correlated hits the requested correlation and moments with empirical = TRUE", {
  set.seed(5)
  d <- sim_correlated(n = 1000, r = 0.3, mu_append = 180, mu_mass = 80,
                      sd_append = 10, sd_mass = 5)
  expect_named(d, c("Append", "Mass"))
  expect_equal(cor(d$Append, d$Mass), 0.3, tolerance = 1e-8)
  expect_equal(mean(d$Append), 180, tolerance = 1e-8)
  expect_equal(sd(d$Mass), 5, tolerance = 1e-8)
})

test_that("transient error in mass alone attenuates the mass-appendage correlation", {
  set.seed(5)
  clean <- sim_correlated(n = 2000, r = 0.3)
  set.seed(5)
  noisy <- sim_correlated(n = 2000, r = 0.3, transient_error_mass = 1)
  expect_lt(abs(cor(noisy$Append, noisy$Mass)), abs(cor(clean$Append, clean$Mass)))
})

test_that("sqrt(calc_lambda) approximates the SMA slope on weakly dispersed traits", {
  set.seed(8)
  d <- sim_allometric(n = 2000, b_avg = 0.33, r_am = 0.3, trim_sd = NULL)
  b_sma <- unname(coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]])
  expect_equal(sqrt(calc_lambda(x = d$Mass, y = d$Append)), b_sma, tolerance = 0.02)
})
