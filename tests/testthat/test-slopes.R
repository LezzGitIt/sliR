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

test_that("averaging slopes that disagree by more than slope_diff_warn warns", {
  d <- make_divergent_slopes()
  expect_warning(build_sli_slopes_tbl(d, control = c("Age", "Sex")), "differ by more than 0.15")

  w <- tryCatch(build_sli_slopes_tbl(d, control = c("Age", "Sex")),
                warning = \(w) conditionMessage(w))
  # The message must name the offending cell and both slopes, not just complain.
  expect_match(w, "Age = ")
  expect_match(w, "b_sma_Age = ")
  expect_match(w, "b_sma_Sex = ")
  expect_match(w, "spread = ")
})

test_that("the warning is silenced by slope_diff_warn = NULL or Inf, and by raising it", {
  d <- make_divergent_slopes()
  expect_no_warning(build_sli_slopes_tbl(d, control = c("Age", "Sex"), slope_diff_warn = NULL))
  expect_no_warning(build_sli_slopes_tbl(d, control = c("Age", "Sex"), slope_diff_warn = Inf))
  expect_no_warning(build_sli_slopes_tbl(d, control = c("Age", "Sex"), slope_diff_warn = 5))
})

test_that("concordant slopes do not warn, and a single control never warns", {
  expect_no_warning(build_sli_slopes_tbl(make_concordant_slopes(), control = c("Age", "Sex")))
  # Nothing is averaged with one control variable, so there is no spread to report.
  expect_no_warning(build_sli_slopes_tbl(make_divergent_slopes(), control = "Age"))
})

test_that("the warning does not alter the returned table", {
  d <- make_divergent_slopes()
  warned <- suppressWarnings(build_sli_slopes_tbl(d, control = c("Age", "Sex")))
  quiet  <- build_sli_slopes_tbl(d, control = c("Age", "Sex"), slope_diff_warn = NULL)
  expect_equal(warned, quiet)
})

test_that("a cell with an NA slope is skipped by the spread check rather than erroring", {
  d <- make_divergent_slopes()
  d$Age[d$Age == "Juv"][1:5] <- NA
  expect_no_error(suppressWarnings(build_sli_slopes_tbl(d, control = c("Age", "Sex"))))
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

## build_sli_slopes_hierarchical() ----

test_that("well-populated cells resolve at the cell level, to their own distinct slope", {
  d   <- dplyr::filter(make_cascade_fixture(), Age != "Rare")
  out <- build_sli_slopes_hierarchical(d, control = c("Age", "Sex"), n_min_cell = 10, n_min_marginal = 15)

  expect_equal(nrow(out), 4L)
  expect_true(all(out$pass_level == "cell"))
  # Four distinct b_sma inputs (0.30/0.40/0.50/0.60) -> four distinct resolved slopes,
  # not a single value that would indicate a fall-through to marginal or pooled.
  expect_equal(length(unique(round(out$b_sli_avg, 2))), 4L)
})

test_that("a cell too sparse for its own fit falls back through marginal to pooled", {
  out <- build_sli_slopes_hierarchical(make_cascade_fixture(), control = c("Age", "Sex"),
                                       n_min_cell = 10, n_min_marginal = 15)

  main <- dplyr::filter(out, Age != "Rare")
  expect_true(all(main$pass_level == "cell"))

  orphan <- dplyr::filter(out, Age == "Rare")
  expect_equal(nrow(orphan), 1L)
  # n = 4 for the cell, the Age marginal, and the Sex marginal alike -- all three fail
  # n_min_cell/n_min_marginal = 10/15, so this is a real pooled fallback, not a coincidence.
  expect_equal(orphan$pass_level, "pooled")

  pooled_slope <- build_sli_slopes_hierarchical(make_cascade_fixture(), control = character(0))$b_sli_avg
  expect_equal(orphan$b_sli_avg, pooled_slope)
})

test_that("an unreliable pooled fit sets every cell's slope to NA, pass_level = 'none'", {
  out <- build_sli_slopes_hierarchical(make_unreliable_fixture(), control = "Sex")
  expect_equal(nrow(out), 2L)
  expect_true(all(out$pass_level == "none"))
  expect_true(all(is.na(out$b_sli_avg)))
})

test_that("control = character(0) returns a single pooled row", {
  d   <- make_cascade_fixture()
  out <- build_sli_slopes_hierarchical(d, control = character(0))
  expect_equal(nrow(out), 1L)
  expect_equal(names(out), c("pass_level", "b_sli_avg"))
  expect_equal(out$pass_level, "pooled")

  fit <- fit_slope_reliability(
    dplyr::mutate(d, .log_app = log(Append), .log_mass = log(Mass)), "Append", "Mass"
  )
  expect_equal(out$b_sli_avg, fit$slope)
})

test_that("a single control variable never resolves at the cell level", {
  out <- build_sli_slopes_hierarchical(make_cascade_fixture(), control = "Age",
                                       n_min_cell = 10, n_min_marginal = 15)
  expect_false("cell" %in% out$pass_level)
  expect_true(all(out$pass_level %in% c("marginal", "pooled")))
})

test_that("more than two control variables errors", {
  d <- make_cascade_fixture()
  d$Extra <- "x"
  expect_error(
    build_sli_slopes_hierarchical(d, control = c("Age", "Sex", "Extra")),
    "length 0, 1, or 2"
  )
})

test_that("NA and every unknown_codes spelling collapse to the same cascade class", {
  d <- make_cascade_fixture()
  rare_idx <- which(d$Age == "Rare" & d$Sex == "X")
  d$Age[rare_idx[1:2]] <- NA
  d$Age[rare_idx[3]]   <- "Unknown"
  # Sex == "X" is entirely NA/unknown-coded Age now; all three spellings must
  # still cascade together as one "Unk" class, not fragment into separate,
  # each-too-small groups.
  out <- build_sli_slopes_hierarchical(d, control = c("Age", "Sex"), n_min_cell = 10, n_min_marginal = 15)

  rare_rows <- dplyr::filter(out, Sex == "X")
  expect_equal(nrow(rare_rows), 3L)  # NA, "Unknown", and the one untouched literal "Rare" row
  expect_equal(length(unique(rare_rows$b_sli_avg)), 1L)
  expect_equal(length(unique(rare_rows$pass_level)), 1L)
})

test_that("calc_sli(method = 'hierarchical') assigns each row its cascade-resolved slope, NAs included", {
  d <- make_cascade_fixture()
  d$Age[1:3] <- NA  # exercise the NA-key join path explicitly

  slopes <- build_sli_slopes_hierarchical(d, control = c("Age", "Sex"), n_min_cell = 10, n_min_marginal = 15)
  out <- calc_sli(d, control = c("Age", "Sex"), method = "hierarchical",
                  n_min_cell = 10, n_min_marginal = 15)

  key <- match(paste(d$Age, d$Sex), paste(slopes$Age, slopes$Sex))
  expected <- d$Append * (mean(d$Mass) / d$Mass)^slopes$b_sli_avg[key]
  expect_equal(out$sli, expected)
  expect_true(anyNA(d$Age))  # confirms the NA-join path was actually exercised, not vacuous
})
