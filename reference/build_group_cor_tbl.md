# Per-group mass–appendage correlation diagnostic

Reports, for each combination of `control` levels, the OLS regression of
mass on appendage length together with its correlation and p-value. Use
it before
[`calc_sli()`](https://LezzGitIt.github.io/sliR/reference/calc_sli.md)`(control = ...)`:
a group whose mass and appendage length are uncorrelated has no
meaningful SMA slope, and should not be used to size-correct its
members.

## Usage

``` r
build_group_cor_tbl(
  df,
  Append = Append,
  Mass = Mass,
  control,
  unknown_codes = c("Unk", "U", "Unknown")
)
```

## Arguments

- df:

  A data frame of individual measurements.

- Append, Mass:

  Unquoted column names for appendage length and body mass. **These
  default to columns literally named `Append` and `Mass`.** If your
  columns are named otherwise, pass them explicitly (e.g.
  `Append = wing, Mass = mass`); an informative error is raised if the
  columns are absent.

- control:

  Character vector of grouping columns. Each must be character or
  factor.

- unknown_codes:

  Values in the `control` columns that mark an unknown group and so
  cannot be assigned a slope under `method = "average"`.

## Value

A tibble with one row per observed combination of `control` levels: the
control columns, `n`, `b_ols` (slope of mass on appendage), `r` (Pearson
correlation), and `p_value`.

## Examples

``` r
set.seed(1)
d <- sim_allometric(n = 400)
d$Sex <- rep(c("F", "M"), length.out = nrow(d))
build_group_cor_tbl(d, control = "Sex")
#> # A tibble: 2 × 5
#>   Sex       n b_ols     r     p_value
#>   <chr> <int> <dbl> <dbl>       <dbl>
#> 1 F       199 0.499 0.353 0.000000318
#> 2 M       198 0.310 0.255 0.000287   
```
