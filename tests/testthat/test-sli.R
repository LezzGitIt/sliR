## A tiny fixture with hand-checkable values: the scalar SLI path is deterministic, so it can be pinned exactly rather than approximately.
toy <- tibble::tibble(
  Append = c(100, 200, 300, 400),
  Mass   = c(10, 20, 30, 40)
)

test_that("calc_sli reproduces the Peig & Green formula exactly", {
  L0  <- mean(toy$Mass) # 25
  out <- calc_sli(toy, b_sli = 0.33)
  expect_equal(out$sli, toy$Append * (L0 / toy$Mass)^0.33)
})

test_that("calc_sli preserves input row order by default and sorts on request", {
  out <- calc_sli(toy, b_sli = 0.33)
  expect_identical(out$Mass, toy$Mass)

  sorted <- calc_sli(toy, b_sli = 0.33, sort = TRUE)
  expect_identical(sorted$sli, sort(sorted$sli, decreasing = TRUE))
})

test_that("M0 rescales SLI by a constant and leaves rankings untouched", {
  default <- calc_sli(toy, b_sli = 0.33)$sli
  custom  <- calc_sli(toy, b_sli = 0.33, M0 = 1000)$sli
  expect_identical(rank(default), rank(custom))
  expect_equal(sd(custom / default), 0)
})

test_that("the reference mass is named M0, not L0", {
  expect_error(calc_sli(toy, b_sli = 0.33, L0 = 1000))
  expect_true("M0" %in% names(formals(calc_sli)))
  expect_false("L0" %in% names(formals(calc_sli)))
})

test_that("b_sli = 0 makes SLI the untransformed appendage length", {
  expect_equal(calc_sli(toy, b_sli = 0)$sli, toy$Append)
})

test_that("rename_col names the output column and NULL leaves it as sli", {
  expect_true("sli" %in% names(calc_sli(toy, b_sli = 0.33)))
  renamed <- calc_sli(toy, b_sli = 0.33, rename_col = "sli_isometry")
  expect_true("sli_isometry" %in% names(renamed))
  expect_false("sli" %in% names(renamed))
})

test_that("calc_sli chains, so several indices can sit side by side", {
  out <- toy |>
    calc_sli(b_sli = 0.33, rename_col = "sli_iso") |>
    calc_sli(b_sli = 0.25, rename_col = "sli_shallow")
  expect_named(out, c("Append", "Mass", "sli_iso", "sli_shallow"))
})

test_that("calc_sli errors informatively when the defaulted columns are absent", {
  lower <- dplyr::rename(toy, wing = Append, mass = Mass)
  expect_error(calc_sli(lower), "not found in `df`")
  expect_no_error(calc_sli(lower, Append = wing, Mass = mass))
})

test_that("calc_sli rejects a malformed rename_col and a non-scalar M0", {
  expect_error(calc_sli(toy, rename_col = c("a", "b")), "single string")
  expect_error(calc_sli(toy, M0 = c(1, 2)), "single number")
})

test_that("calc_sli warns that b_sli is ignored once control is supplied", {
  d <- sim_allometric(n = 100)
  d$Sex <- rep(c("F", "M"), length.out = nrow(d))
  expect_warning(calc_sli(d, b_sli = 0.4, control = "Sex"), "ignored")
})

test_that("calc_sli forwards slope_diff_warn to build_sli_slopes_tbl", {
  d <- make_divergent_slopes()
  expect_warning(calc_sli(d, control = c("Age", "Sex")), "differ by more than")
  expect_no_warning(calc_sli(d, control = c("Age", "Sex"), slope_diff_warn = NULL))
})
