## A grouped fixture: sexes differ in their allometric slope, so per-group SMA slopes should differ too.
make_grouped <- function(n = 400, seed = 11) {
  set.seed(seed)
  d <- sim_allometric(n = n)
  d$Sex <- rep(c("F", "M"), length.out = nrow(d))
  d$Age <- factor(rep(c("Juv", "Adult"), each = 2, length.out = nrow(d)))
  d
}

test_that("build_sli_slopes_tbl returns one row per observed level combination", {
  d   <- make_grouped()
  out <- build_sli_slopes_tbl(d, control = c("Age", "Sex"))
  expect_equal(nrow(out), 4L)
  expect_true(all(c("Age", "Sex", "n", "b_sma_Age", "b_sma_Sex", "b_sli_avg") %in% names(out)))
  expect_equal(sum(out$n), nrow(d))
})

test_that("b_sli_avg is the row mean of the per-variable SMA slopes", {
  out <- build_sli_slopes_tbl(make_grouped(), control = c("Age", "Sex"))
  expect_equal(out$b_sli_avg, rowMeans(cbind(out$b_sma_Age, out$b_sma_Sex)))
})

test_that("build_sli_slopes_tbl preserves the key column's type for the join back", {
  out <- build_sli_slopes_tbl(make_grouped(), control = c("Age", "Sex"))
  expect_s3_class(out$Age, "factor")
  expect_type(out$Sex, "character")
})

test_that("unknown-coded and NA groups are dropped from the slope table", {
  d <- make_grouped()
  d$Sex[1:10] <- "Unk"
  out <- build_sli_slopes_tbl(d, control = "Sex")
  expect_false("Unk" %in% out$Sex)
  expect_equal(sum(out$n), nrow(d) - 10)
})

test_that("build_sli_slopes_tbl rejects a numeric grouping column", {
  d <- make_grouped()
  d$grp <- seq_len(nrow(d))
  expect_error(build_sli_slopes_tbl(d, control = "grp"), "character or factor")
})

test_that("calc_sli(control=) assigns each row its group's averaged slope", {
  d      <- make_grouped()
  slopes <- build_sli_slopes_tbl(d, control = "Sex")
  out    <- calc_sli(d, control = "Sex")

  L0       <- mean(d$Mass)
  expected <- d$Append * (L0 / d$Mass)^slopes$b_sli_avg[match(d$Sex, slopes$Sex)]
  expect_equal(out$sli, expected)
})

test_that("calc_sli(control=) gives NA to unknown and missing groups, and keeps every row", {
  d <- make_grouped()
  d$Sex[1:10] <- "Unk"
  d$Age[11:15] <- NA

  out <- calc_sli(d, control = c("Age", "Sex"))
  expect_equal(nrow(out), nrow(d))
  expect_identical(out$Mass, d$Mass)
  expect_equal(sum(is.na(out$sli)), 15L)
  expect_false("b_sli_avg" %in% names(out))
})

test_that("build_group_cor_tbl reports n, slope, correlation, and p per group", {
  out <- build_group_cor_tbl(make_grouped(), control = "Sex")
  expect_equal(nrow(out), 2L)
  expect_named(out, c("Sex", "n", "b_ols", "r", "p_value"))
  expect_false(dplyr::is_grouped_df(out))
  expect_true(all(out$r >= -1 & out$r <= 1))
  expect_true(all(out$p_value >= 0 & out$p_value <= 1))
})

test_that("build_group_cor_tbl agrees with a hand-fitted lm on one group", {
  d   <- make_grouped()
  out <- build_group_cor_tbl(d, control = "Sex")
  f   <- dplyr::filter(d, Sex == "F")
  mod <- lm(Mass ~ Append, data = f)

  row <- dplyr::filter(out, Sex == "F")
  expect_equal(row$b_ols, unname(coef(mod)[2]))
  expect_equal(row$r, cor(f$Append, f$Mass))
})

test_that("build_group_cor_tbl works when trait columns are lowercase", {
  d <- dplyr::rename(make_grouped(), wing = Append, mass = Mass)
  expect_no_error(build_group_cor_tbl(d, Append = wing, Mass = mass, control = "Sex"))
})
