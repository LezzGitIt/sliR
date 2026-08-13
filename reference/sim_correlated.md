# Simulate correlated appendage and mass on the raw scale

Draws a bivariate normal pair of appendage length and body mass with a
given correlation, then optionally perturbs each with measurement error
and transient biological fluctuation. Unlike
[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)
there is no log transform, no allometric target, and no gradient: this
is the tool for building intuition about how a fitted slope responds to
correlation and to error added on one trait but not the other.

## Usage

``` r
sim_correlated(
  n = 3000,
  r_app_mass = 0.3,
  mu_append = 180,
  mu_mass = 80,
  sd_append = 10,
  sd_mass = 5,
  meas_error = 0,
  transient_error_append = 0,
  transient_error_mass = 0,
  empirical = TRUE,
  r = NULL
)
```

## Arguments

- n:

  Number of individuals.

- r_app_mass:

  Correlation between appendage length and mass.

- mu_append, mu_mass:

  Means of the two traits.

- sd_append, sd_mass:

  Standard deviations of the two traits.

- meas_error:

  Observer error, as a fraction of each trait's standard deviation.
  Applied to both traits.

- transient_error_append, transient_error_mass:

  Biological fluctuation, as a fraction of that trait's standard
  deviation.

- empirical:

  If `TRUE` (the default), the drawn sample has exactly the requested
  moments and correlations, rather than being a random draw from a
  population with them. Trimming and error break this guarantee.

- r:

  Soft-deprecated: use `r_app_mass` instead, the name
  [`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)
  uses for the same quantity.

## Value

A tibble of `n` rows with columns `Append` and `Mass`.

## Examples

``` r
set.seed(1)
# Transient fluctuation in mass alone attenuates the fitted slope
clean <- sim_correlated(r_app_mass = 0.3, transient_error_mass = 0)
noisy <- sim_correlated(r_app_mass = 0.3, transient_error_mass = 1)
c(clean = cor(clean$Append, clean$Mass), noisy = cor(noisy$Append, noisy$Mass))
#>     clean     noisy 
#> 0.3000000 0.1978644 
```
