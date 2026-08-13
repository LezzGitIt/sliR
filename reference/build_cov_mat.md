# Target covariance matrix for appendage, mass, and an optional gradient

Assembles the log-scale covariance matrix implied by an allometry,
optionally extended with an environmental gradient. The allometry is
resolved by
[`implied_allometry()`](https://LezzGitIt.github.io/sliR/reference/implied_allometry.md),
so supply any three of `r_app_mass`, `b_ols`, `b_sma`, `b_avg`,
`sd_log_append`, `sd_log_mass`.

## Usage

``` r
build_cov_mat(
  r_app_mass = NULL,
  b_ols = NULL,
  b_sma = NULL,
  b_avg = NULL,
  sd_log_append = NULL,
  sd_log_mass = NULL,
  gradient = NULL,
  r_grad_app = 0,
  r_grad_mass = 0,
  sd_gradient = 1
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

- gradient:

  Optional string naming an environmental gradient (e.g.
  `"Temperature"`). `NULL` (the default) returns the 2x2 morphological
  block.

- r_grad_app, r_grad_mass:

  Correlations of the gradient with `log(Append)` and `log(Mass)`.
  Ignored when `gradient` is `NULL`.

- sd_gradient:

  Standard deviation of the gradient. Defaults to `1`. Ignored when
  `gradient` is `NULL`. The morphological correlations and slopes do not
  depend on it; it scales only the gradient's own variance and
  covariances.

## Value

A covariance matrix, 2x2 with dimnames `Append`, `Mass`, or 3x3 with the
gradient's name appended.

## Details

The gradient enters at **unit standard deviation by default**, because
this matrix primarily describes a correlation structure and
[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)
rescales the drawn gradient itself. Set `sd_gradient` to put the
gradient on its own scale — e.g. `sd_gradient = 0.18` for a temperature
index — so the returned matrix is the covariance of your actual system
rather than a standardised one.

## See also

[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md),
which draws from this structure.

## Examples

``` r
build_cov_mat(b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07)
#>              Append    Mass
#> Append 0.0005444444 0.00049
#> Mass   0.0004900000 0.00490

# A temperature gradient on its own scale (SD 0.18), not standardised
build_cov_mat(
  b_sma = 1/3, r_app_mass = 0.3, sd_log_mass = 0.07,
  gradient = "Temperature", r_grad_app = -0.3, r_grad_mass = -0.1,
  sd_gradient = 0.18
)
#>                    Append     Mass Temperature
#> Append       0.0005444444  0.00049    -0.00126
#> Mass         0.0004900000  0.00490    -0.00126
#> Temperature -0.0012600000 -0.00126     0.03240
```
