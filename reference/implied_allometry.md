# Resolve an allometry from any sufficient subset of its parameters

A log-log allometry between an appendage and body mass has three free
parameters: the two log-scale standard deviations and the correlation
between them. Six quantities describe it, and they are linked by

## Usage

``` r
implied_allometry(
  r_app_mass = NULL,
  b_ols = NULL,
  b_sma = NULL,
  b_avg = NULL,
  sd_log_append = NULL,
  sd_log_mass = NULL
)
```

## Arguments

- r_app_mass:

  Pearson correlation between `log(Append)` and `log(Mass)`. Must be
  non-zero.

- b_ols, b_sma, b_avg:

  The OLS, SMA, and midpoint slopes of `log(Append)` on `log(Mass)`.

- sd_log_append, sd_log_mass:

  Log-scale standard deviations.

## Value

A one-row tibble with all six quantities.

## Details

\$\$b\_{SMA} = \mathrm{sign}(r)\\\sigma_A / \sigma_M, \quad b\_{OLS} = r
\\ b\_{SMA}, \quad b\_{avg} = (b\_{OLS} + b\_{SMA}) / 2 .\$\$

So *any two* of `r_app_mass`, `b_ols`, `b_sma`, `b_avg` determine the
other two, and one standard deviation then fixes the scale.

Supply **two shape quantities, or none**. Supplying none takes the
isometric default, `b_sma = 1/3, r_app_mass = 0.3`. Supplying exactly
one is an error: pairing your slope with a default correlation would
silently determine the other two slopes. Over-specification is allowed
but checked for consistency.

The standard deviations behave differently. Neither slope nor
correlation depends on them, so they take an ordinary default
(`sd_log_mass = 0.07`) when omitted. They set raw-scale dispersion and
skew, which still matters, because the SLI is computed on raw metrics.

The three slopes are not interchangeable, and the difference matters:

- `b_ols` is the slope of \\E\[\log A \mid \log M\]\\. In the bivariate
  lognormal that
  [`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)
  draws from, this is the *generative* slope: the coefficient an
  [`lm()`](https://rdrr.io/r/stats/lm.html) will recover.

- `b_sma` is \\\sigma_A / \sigma_M\\, the slope a standardised major
  axis fit returns. It is the exponent
  [`calc_sli()`](https://LezzGitIt.github.io/sliR/reference/calc_sli.md)
  assumes and
  [`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md)
  estimates, and the one usually meant by "the allometric slope".

- `b_avg` is their midpoint. No estimator returns it. It is a convention
  for splitting the difference when OLS and SMA disagree, and is
  retained because earlier work used it. **Note that `b_avg = 1/3` does
  not mean isometry**: with `b_avg = 0.33, r_app_mass = 0.3` neither
  slope is near 1/3.

## Examples

``` r
# Isometry as an SMA fit would report it
implied_allometry(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07)
#> # A tibble: 1 × 6
#>   r_app_mass b_ols b_sma b_avg sd_log_append sd_log_mass
#>        <dbl> <dbl> <dbl> <dbl>         <dbl>       <dbl>
#> 1        0.3   0.1 0.333 0.217        0.0233        0.07

# The parameterisation used by the SLI manuscript
implied_allometry(b_avg = 0.33, r_app_mass = 0.3, sd_log_mass = 0.07)
#> # A tibble: 1 × 6
#>   r_app_mass b_ols b_sma b_avg sd_log_append sd_log_mass
#>        <dbl> <dbl> <dbl> <dbl>         <dbl>       <dbl>
#> 1        0.3 0.152 0.508  0.33        0.0355        0.07

# Pin both slopes and let the correlation follow
implied_allometry(b_ols = 0.15, b_sma = 0.50, sd_log_mass = 0.07)
#> # A tibble: 1 × 6
#>   r_app_mass b_ols b_sma b_avg sd_log_append sd_log_mass
#>        <dbl> <dbl> <dbl> <dbl>         <dbl>       <dbl>
#> 1        0.3  0.15   0.5 0.325         0.035        0.07
```
