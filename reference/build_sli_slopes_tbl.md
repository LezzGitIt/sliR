# Per-group SMA allometric slopes

Fits one standardised major axis (SMA) regression of `log(Append)` on
`log(Mass)` per `control` variable, letting the slope vary across that
variable's levels. Every combination of levels present in the data then
receives `b_sli_avg`, the mean of the slopes its levels attract.

## Usage

``` r
build_sli_slopes_tbl(
  df,
  Append = Append,
  Mass = Mass,
  control,
  unknown_codes = c("Unk", "U", "Unknown"),
  slope_diff_warn = 0.15
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

- slope_diff_warn:

  Warn when the per-variable SMA slopes being averaged within a cell
  span more than this. Defaults to `0.15`, roughly half the isometric
  exponent. `NULL` or `Inf` disables the check. Has no effect with a
  single `control` variable, where nothing is averaged.

## Value

A tibble with one row per observed combination of `control` levels: the
control columns, `n` (rows contributing to that cell), one `b_sma_<var>`
column per control variable, and their row mean `b_sli_avg`. A cell
whose slope is missing for any control variable gets `b_sli_avg = NA`.

## Details

Note the structure: with `control = c("Age", "Sex")` this fits **two**
models — one interacting mass with age, one interacting mass with sex —
and averages their slopes per age × sex cell. It does not fit a single
age × sex interaction (which is not possible in the smatr package). That
keeps the per-cell sample sizes at the marginal rather than the joint
level, which matters when some cells are sparse.

The cost of that averaging is that a cell whose control variables
disagree about the slope is collapsed to their midpoint, and nothing in
`b_sli_avg` records the disagreement. `slope_diff_warn` guards against
this: when any cell's per-variable slopes span more than that much, a
warning names the worst cell and its spread. Inspect the `b_sma_<var>`
columns before trusting `b_sli_avg`, and consider correcting on one
control variable at a time.

Rows whose `control` values are `NA` or match `unknown_codes` are
excluded from slope fitting and from the returned table.

## See also

[`build_group_cor_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_group_cor_tbl.md),
which flags groups whose mass–appendage relationship is too weak for a
per-group slope to be meaningful.

## Examples

``` r
set.seed(1)
d <- sim_allometric(n = 400)
d$Sex <- rep(c("F", "M"), length.out = nrow(d))
build_sli_slopes_tbl(d, control = "Sex")
#> # A tibble: 2 × 4
#>   Sex       n b_sma_Sex b_sli_avg
#>   <chr> <int>     <dbl>     <dbl>
#> 1 F       199     0.314     0.314
#> 2 M       198     0.366     0.366
```
