# Simulate log-normal morphology, optionally along an environmental gradient

Draws appendage length and body mass from a bivariate log-normal with a
specified allometry, and optionally an exogenous environmental gradient
(temperature, latitude, year, rainfall) that both traits respond to
linearly.

## Usage

``` r
sim_allometric(
  n = 3000,
  r_app_mass = NULL,
  b_ols = NULL,
  b_sma = NULL,
  b_avg = NULL,
  sd_log_append = NULL,
  sd_log_mass = NULL,
  mean_append = 180,
  mean_mass = 80,
  gradient = NULL,
  gradient_dist = c("uniform", "normal"),
  gradient_range = NULL,
  mean_gradient = 0,
  sd_gradient = 1,
  r_grad_app = 0,
  r_grad_mass = 0,
  meas_error = 0,
  transient_error_append = 0,
  transient_error_mass = 0,
  trim_sd = 3,
  log = TRUE,
  empirical = TRUE
)
```

## Arguments

- n:

  Number of individuals to draw, before trimming.

- r_app_mass:

  Pearson correlation between `log(Append)` and `log(Mass)`. Must be
  non-zero.

- b_ols, b_sma, b_avg:

  The OLS, SMA, and midpoint slopes of `log(Append)` on `log(Mass)`.

- sd_log_append, sd_log_mass:

  Log-scale standard deviations.

- mean_append, mean_mass:

  Raw-scale means, used as the log-scale location via
  [`log()`](https://rdrr.io/r/base/Log.html).

- gradient:

  Optional string naming an environmental gradient (e.g.
  `"Temperature"`). `NULL` (the default) returns the 2x2 morphological
  block.

- gradient_dist:

  Marginal distribution of the gradient: `"uniform"` (the default) or
  `"normal"`.

- gradient_range:

  Two increasing numbers giving the gradient's bounds, e.g.
  `c(1970, 2026)`. Uniform gradients only. Determines `mean_gradient`
  and `sd_gradient`, which must then not be supplied.

- mean_gradient, sd_gradient:

  Location and spread of the gradient. Default to `0` and `1`, a
  z-scored gradient.

- r_grad_app, r_grad_mass:

  Correlations of the gradient with `log(Append)` and `log(Mass)`.
  Ignored when `gradient` is `NULL`.

- meas_error:

  Observer error, as a fraction of each trait's standard deviation.
  Applied to both traits.

- transient_error_append, transient_error_mass:

  Biological fluctuation, as a fraction of that trait's standard
  deviation.

- trim_sd:

  Drop individuals more than `trim_sd` raw-scale standard deviations
  from the mean of `Append` or `Mass`. `NULL` disables trimming. The
  gradient is never trimmed. Because the traits are log-normal,
  raw-scale trimming is asymmetric and removes more of the right tail
  than the left.

- log:

  If `TRUE` (the default), append the `Append_log` and `Mass_log`
  columns, computed after error is added. `FALSE` returns the raw traits
  and the gradient only.

- empirical:

  If `TRUE` (the default), the drawn sample has exactly the requested
  moments and correlations, rather than being a random draw from a
  population with them. Trimming and error break this guarantee.

## Value

A tibble with `n` rows or fewer: `Append`, `Mass`, the gradient column
if requested, and `Append_log`, `Mass_log` unless `log = FALSE`.

## Details

The allometry is resolved by
[`implied_allometry()`](https://LezzGitIt.github.io/sliR/reference/implied_allometry.md).
Supply **two** of `r_app_mass`, `b_ols`, `b_sma`, `b_avg`, or none to
accept the isometric default `b_sma = 1/3, r_app_mass = 0.3`. Supplying
exactly one is an error. `sd_log_append` and `sd_log_mass` set the scale
and default to `sd_log_mass = 0.07`; they leave every slope and
correlation unchanged.

## The gradient

The gradient is a **sampling design variable**, not a third trait, so it
is drawn exogenously and the log traits are built linearly from it. Its
marginal is therefore free, and it defaults to `"uniform"` because:

- **it is bounded**, so `gradient_range = c(-90, 90)` really does yield
  latitudes within \\\pm 90\\, and years stay inside their range;

- **it covers the gradient evenly.** Across 15 equal-width bins a normal
  gradient leaves the extreme bins nearly empty, with the fullest bin
  thousands of times the emptiest; a uniform one is flat.

`gradient_dist = "normal"` recovers the ordinary three-variable
multivariate normal. Either way `r_grad_app` is the correlation of the
gradient with `log(Append)` (not with `Append`), and \\E\[\log A \mid
G\]\\ is linear in \\G\\ by construction.

Not every set of correlations is attainable. An appendage that lengthens
while mass falls along the gradient is incompatible with a strong
positive appendage-mass correlation, and such a request is rejected
rather than silently approximated.

## Error

Error is added on the raw scale, after exponentiation, each component a
fraction of the trait's own standard deviation. `meas_error` is observer
error and applies to both traits alike; `transient_error_*` represents
real short-term biological fluctuation, which afflicts mass far more
than a skeletal appendage.

## Reproducibility and trimming

This function consumes the random number stream, so call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before it for
reproducible output. Because the internal order of the draws is an
implementation detail, byte-for-byte reproducibility is tied to a fixed
package version: pin one (e.g.
`remotes::install_github("LezzGitIt/sliR@v0.1.0")`) for analyses that
must reproduce exactly.

Trimming (`trim_sd`) is applied **jointly** to `Append` and `Mass`: a
row is dropped if either trait exceeds the cut. Code that instead trims
one trait and then the other, using each trait's pre-trim standard
deviation, will retain a slightly different set of rows for the same
draw.

## See also

[`implied_allometry()`](https://LezzGitIt.github.io/sliR/reference/implied_allometry.md)
to inspect what a parameter set implies, and
[`sim_correlated()`](https://LezzGitIt.github.io/sliR/reference/sim_correlated.md)
for a raw-scale bivariate draw with no allometric target.

## Examples

``` r
set.seed(1)

# Isometry, no gradient
d <- sim_allometric(n = 500)
coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]]
#> [1] 0.3280054

# Allen's rule: warmer places, relatively longer appendages
warm <- sim_allometric(
  n = 500, gradient = "Temperature", gradient_range = c(0, 24),
  r_grad_app = 0.3, r_grad_mass = -0.1
)
range(warm$Temperature)
#> [1]  0.03595807 23.96404193

# A century of specimens, drawn evenly across years
yr <- sim_allometric(n = 500, gradient = "Year", gradient_range = c(1970, 2026))

# Raw traits only
names(sim_allometric(n = 10, log = FALSE))
#> [1] "Append" "Mass"  
```
