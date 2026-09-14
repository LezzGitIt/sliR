## The reference allometry every combination below must reproduce.
target <- list(r_app_mass = 0.3, b_ols = 0.1523077, b_sma = 0.5076923, b_avg = 0.33)

test_that("any two shape quantities plus one SD recover the same allometry", {
  pairs <- list(
    c("r_app_mass", "b_sma"), c("r_app_mass", "b_ols"), c("r_app_mass", "b_avg"),
    c("b_ols", "b_sma"),      c("b_ols", "b_avg"),      c("b_sma", "b_avg")
  )
  for (pr in pairs) {
    out <- do.call(implied_allometry, c(target[pr], list(sd_log_mass = 0.07)))
    expect_equal(out$r_app_mass, target$r_app_mass, tolerance = 1e-5, info = paste(pr, collapse = "+"))
    expect_equal(out$b_ols,      target$b_ols,      tolerance = 1e-5, info = paste(pr, collapse = "+"))
    expect_equal(out$b_sma,      target$b_sma,      tolerance = 1e-5, info = paste(pr, collapse = "+"))
    expect_equal(out$b_avg,      target$b_avg,      tolerance = 1e-5, info = paste(pr, collapse = "+"))
  }
})

test_that("the scale can be pinned from either standard deviation", {
  by_mass <- implied_allometry(b_sma = 0.5, r_app_mass = 0.3, sd_log_mass = 0.07)
  expect_equal(by_mass$sd_log_append, 0.5 * 0.07)

  by_app <- implied_allometry(b_sma = 0.5, r_app_mass = 0.3, sd_log_append = 0.035)
  expect_equal(by_app$sd_log_mass, 0.035 / 0.5)
})

test_that("one shape quantity plus both SDs is sufficient", {
  out <- implied_allometry(r_app_mass = 0.3, sd_log_append = 0.035, sd_log_mass = 0.07)
  expect_equal(out$b_sma, 0.5)
  expect_equal(out$b_ols, 0.15)
})

test_that("b_sma alongside both SDs is under-determined, since the SDs already say b_sma", {
  expect_error(
    implied_allometry(b_sma = 0.5, sd_log_append = 0.035, sd_log_mass = 0.07),
    "under-determined"
  )
})

test_that("exactly one shape quantity errors rather than pairing with a default correlation", {
  expect_error(implied_allometry(b_sma = 1/3), "under-determined")
  expect_error(implied_allometry(b_ols = 0.15), "under-determined")
  expect_error(implied_allometry(b_avg = 0.33, sd_log_mass = 0.07), "under-determined")
})

test_that("two shape quantities suffice; the scale takes an ordinary default", {
  out <- implied_allometry(b_sma = 1/3, r_app_mass = 0.3)
  expect_equal(out$sd_log_mass, 0.07)
  expect_equal(out$b_sma, 1/3)
})

test_that("supplying nothing uses the documented isometric default", {
  out <- implied_allometry()
  expect_equal(out$b_sma, 1/3)
  expect_equal(out$r_app_mass, 0.3)
  expect_equal(out$sd_log_mass, 0.07)
  expect_equal(out$b_ols, 0.1)
})

test_that("a scale-only call keeps the isometric default shape", {
  out <- implied_allometry(sd_log_mass = 0.10)
  expect_equal(out$b_sma, 1/3)
  expect_equal(out$sd_log_mass, 0.10)
})

test_that("scale changes no slope and no correlation", {
  a <- implied_allometry(b_avg = 0.33, r_app_mass = 0.3, sd_log_mass = 0.05)
  b <- implied_allometry(b_avg = 0.33, r_app_mass = 0.3, sd_log_mass = 0.30)
  expect_equal(a$b_ols, b$b_ols)
  expect_equal(a$b_sma, b$b_sma)
  expect_equal(a$r_app_mass, b$r_app_mass)
})

test_that("consistent over-specification passes; inconsistent over-specification errors", {
  expect_no_error(
    implied_allometry(b_ols = 0.1523077, b_sma = 0.5076923, r_app_mass = 0.3, sd_log_mass = 0.07)
  )
  expect_error(
    implied_allometry(b_ols = 0.2, b_sma = 0.5, r_app_mass = 0.9, sd_log_mass = 0.07),
    "over-determined and inconsistent"
  )
})

test_that("infeasible parameter sets are rejected", {
  # |b_ols| > |b_sma| implies |r| > 1
  expect_error(implied_allometry(b_ols = 0.6, b_sma = 0.5, sd_log_mass = 0.07), "outside \\[-1, 1\\]")
  expect_error(implied_allometry(r_app_mass = 0, b_sma = 1/3, sd_log_mass = 0.07), "non-zero")
  expect_error(implied_allometry(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = -1), "must be positive")
})

test_that("solver arguments must be single numbers", {
  expect_error(implied_allometry(b_sma = c(0.3, 0.4), r_app_mass = 0.3, sd_log_mass = 0.07), "single non-missing number")
  expect_error(implied_allometry(b_sma = NA, r_app_mass = 0.3, sd_log_mass = 0.07), "single non-missing number")
})

test_that("the manuscript's parameterisation still resolves to its known values", {
  out <- implied_allometry(b_avg = 0.33, r_app_mass = 0.3, sd_log_mass = 0.07)
  expect_equal(out$b_ols, 0.152308, tolerance = 1e-5)
  expect_equal(out$b_sma, 0.507692, tolerance = 1e-5)
})

test_that("b_avg = 1/3 is not isometry, and b_sma = 1/3 is", {
  # Under b_avg = 0.33 neither fitted slope is near 1/3 -- the point of the b_sma default.
  mid <- implied_allometry(b_avg = 0.33, r_app_mass = 0.3, sd_log_mass = 0.07)
  expect_false(isTRUE(all.equal(mid$b_sma, 1/3, tolerance = 0.05)))
  expect_false(isTRUE(all.equal(mid$b_ols, 1/3, tolerance = 0.05)))

  iso <- implied_allometry(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07)
  expect_equal(iso$b_sma, 1/3)
})

test_that("implied_gradient_effect matches the SLI-allens-rule manuscript's beta_T formula", {
  # b_avg = 0.33, r_app_mass = 0.3 implies b_sma = 0.5076923 (shared fixture, test-solve.R:2).
  out <- implied_gradient_effect(b_avg = 0.33, r_app_mass = 0.3,
                                 r_grad_app = -0.3, r_grad_mass = -0.5)
  b_sma <- 0.5076923
  expect_equal(out$b_sma, b_sma, tolerance = 1e-5)
  # Manuscript: beta_T = r_13 * b_sma - 0.33 * r_23 (r_13 = r_grad_app, r_23 = r_grad_mass here).
  expect_equal(out$beta_ref, -0.3 * b_sma - (1 / 3) * -0.5, tolerance = 1e-5)
})

test_that("the two-pathways decomposition sums exactly to beta_ref", {
  out <- implied_gradient_effect(b_sma = 0.5, r_app_mass = 0.3,
                                 r_grad_app = 0.2, r_grad_mass = -0.1, b_anchor = 0.4)
  expect_equal(out$allometry_component + out$differential_component, out$beta_ref)
})

test_that("anchoring at the realised allometry zeroes allometry_component, not differential_component", {
  # Mirrors calc_sli(control = ...): the exponent IS the group's own b_sma.
  out <- implied_gradient_effect(b_sma = 0.5, r_app_mass = 0.3,
                                 r_grad_app = 0.2, r_grad_mass = -0.1, b_anchor = 0.5)
  expect_equal(out$allometry_component, 0)
  expect_false(isTRUE(all.equal(out$differential_component, 0)))
  expect_equal(out$beta_ref, out$differential_component)
})

test_that("no differential association leaves only allometry_component", {
  out <- implied_gradient_effect(b_sma = 0.5, r_app_mass = 0.3,
                                 r_grad_app = 0.2, r_grad_mass = 0.2)
  expect_equal(out$differential_component, 0)
  expect_equal(out$beta_ref, out$allometry_component)
})

test_that("b_anchor defaults to isometry, matching implied_allometry's own default shape", {
  out <- implied_gradient_effect(r_grad_app = 0.1, r_grad_mass = -0.1)  # no shape args -> isometric default
  expect_equal(out$b_anchor, 1 / 3)
  expect_equal(out$b_sma, 1 / 3)  # isometric default shape, as implied_allometry() documents
})

test_that("implied_gradient_effect rejects an infeasible three-way correlation triple", {
  expect_error(
    implied_gradient_effect(b_sma = 1 / 3, r_app_mass = 0.9,
                            r_grad_app = 0.9, r_grad_mass = -0.9),
    "not positive definite"
  )
})

test_that("implied_gradient_effect inherits implied_allometry's own validation", {
  expect_error(
    implied_gradient_effect(b_sma = 0.5, r_grad_app = 0.1, r_grad_mass = -0.1),
    "under-determined"
  )
})

test_that("gradient correlations and the anchor must be single numbers", {
  expect_error(
    implied_gradient_effect(b_sma = 1 / 3, r_app_mass = 0.3, r_grad_app = c(0.1, 0.2), r_grad_mass = -0.1),
    "single non-missing number"
  )
  expect_error(
    implied_gradient_effect(b_sma = 1 / 3, r_app_mass = 0.3, r_grad_app = 0.1, r_grad_mass = -0.1, b_anchor = NA),
    "single non-missing number"
  )
})

test_that("sim_allometric realises whichever slope was pinned", {
  set.seed(1)
  d <- sim_allometric(n = 3000, b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07, trim_sd = NULL)
  expect_equal(unname(coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]]), 1/3, tolerance = 1e-6)

  set.seed(1)
  d2 <- sim_allometric(n = 3000, b_ols = 0.15, r_app_mass = 0.3, sd_log_mass = 0.07, trim_sd = NULL)
  expect_equal(unname(coef(stats::lm(Append_log ~ Mass_log, data = d2))[2]), 0.15, tolerance = 1e-6)
})
