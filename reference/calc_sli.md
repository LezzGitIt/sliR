# Standardized Length Index (SLI)

Size-corrects an appendage length by rescaling every individual to a
common reference body mass, following the scaled mass index of Peig &
Green (2009). Whereas their index standardises *mass* to a reference
*length*, the SLI standardises a *length* to a reference *mass*:

## Usage

``` r
calc_sli(
  df,
  Append = Append,
  Mass = Mass,
  b_sli = 0.33,
  M0 = NULL,
  control = NULL,
  unknown_codes = c("Unk", "U", "Unknown"),
  slope_diff_warn = 0.15,
  rename_col = NULL,
  sort = FALSE
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

- b_sli:

  Scalar allometric scaling exponent. The default `0.33` (\\\approx
  1/3\\) assumes isometry, under which linear dimensions scale with the
  cube root of mass. Ignored when `control` is supplied. Set it to an
  empirically estimated SMA slope to relax isometry.

- M0:

  Reference mass to which all individuals are standardised. `NULL` (the
  default) uses the mean of `Mass`. \\M_0\\ is an arbitrary scalar: it
  rescales every SLI value by a constant factor and so changes the units
  of the index, never the ranking of individuals or the result of any
  subsequent location-invariant analysis. (Peig & Green call this
  \\L_0\\, because their index standardises to a reference *length*; the
  SLI inverts that.)

- control:

  Optional character vector of grouping columns (e.g.
  `c("Age", "Sex")`). When supplied,
  [`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md)
  estimates a separate SMA slope per group and each individual receives
  their group's averaged slope in place of `b_sli`. Individuals whose
  `control` values are `NA` or match `unknown_codes` receive `sli = NA`.

- unknown_codes:

  Values in the `control` columns that mark an unknown group and so
  cannot be assigned a slope.

- slope_diff_warn:

  Passed to
  [`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md).
  Warn when the per-variable SMA slopes averaged within a group differ
  by more than this. Only relevant with two or more `control` variables.

- rename_col:

  Optional string naming the output column. `NULL` (the default) leaves
  it as `sli`.

- sort:

  If `TRUE`, return rows sorted by descending SLI. Defaults to `FALSE`,
  preserving the input row order.

## Value

`df` with one additional numeric column (`sli`, or `rename_col`).

## Details

\$\$SLI_i = A_i \left( \frac{M_0}{M_i} \right)^{b}\$\$

where \\A_i\\ and \\M_i\\ are individual \\i\\'s appendage length and
body mass, \\M_0\\ is an arbitrary reference mass, and \\b\\ is the
allometric scaling exponent relating length to mass.

## References

Peig, J. & Green, A.J. (2009) New perspectives for estimating body
condition from mass/length data: the scaled mass index as an alternative
method. *Oikos* 118, 1883–1891.
[doi:10.1111/j.1600-0706.2009.17643.x](https://doi.org/10.1111/j.1600-0706.2009.17643.x)

## See also

[`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md)
for the per-group slopes, and
[`build_group_cor_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_group_cor_tbl.md)
to check whether per-group slopes are warranted.

## Examples

``` r
set.seed(1)
d <- sim_allometric(n = 200)

# Isometric standardisation
head(calc_sli(d, b_sli = 0.33))
#> # A tibble: 6 × 5
#>   Append  Mass Append_log Mass_log   sli
#>    <dbl> <dbl>      <dbl>    <dbl> <dbl>
#> 1   177.  77.7       5.18     4.35  179.
#> 2   179.  71.5       5.19     4.27  186.
#> 3   175.  71.5       5.16     4.27  182.
#> 4   187.  82.8       5.23     4.42  185.
#> 5   183.  94.1       5.21     4.54  174.
#> 6   174.  67.2       5.16     4.21  185.

# Two indices side by side, in the original row order
d |>
  calc_sli(b_sli = 0.33, rename_col = "sli_isometry") |>
  calc_sli(b_sli = 0.25, rename_col = "sli_shallow") |>
  head()
#> # A tibble: 6 × 6
#>   Append  Mass Append_log Mass_log sli_isometry sli_shallow
#>    <dbl> <dbl>      <dbl>    <dbl>        <dbl>       <dbl>
#> 1   177.  77.7       5.18     4.35         179.        178.
#> 2   179.  71.5       5.19     4.27         186.        184.
#> 3   175.  71.5       5.16     4.27         182.        180.
#> 4   187.  82.8       5.23     4.42         185.        186.
#> 5   183.  94.1       5.21     4.54         174.        176.
#> 6   174.  67.2       5.16     4.21         185.        182.
```
