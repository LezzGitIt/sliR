# Resolve the ground-truth effect of a gradient on relative appendage size

Extends
[`implied_allometry()`](https://LezzGitIt.github.io/sliR/reference/implied_allometry.md)
to a third variable – an environmental or temporal gradient along which
both appendage length and mass may covary (temperature, latitude,
year...). Answers the question a simulation needs before it can score
any estimator: *given this allometry and these gradient correlations,
what relative change should a correctly-specified estimator recover?*

## Usage

``` r
implied_gradient_effect(
  r_app_mass = NULL,
  b_ols = NULL,
  b_sma = NULL,
  b_avg = NULL,
  r_grad_app,
  r_grad_mass,
  b_anchor = 1/3,
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

- r_grad_app:

  Pearson correlation between the gradient and `log(Append)`.

- r_grad_mass:

  Pearson correlation between the gradient and `log(Mass)`.

- b_anchor:

  The scaling exponent defining "no relative change" – the reference
  [`calc_sli()`](https://LezzGitIt.github.io/sliR/reference/calc_sli.md)-style
  methods are implicitly or explicitly scored against. Defaults to
  `1/3`, the isometric exponent for a linear appendage against a
  volumetric body-size proxy like mass.

- sd_log_append, sd_log_mass:

  Log-scale standard deviations.

## Value

A one-row tibble: the resolved allometry (as
[`implied_allometry()`](https://LezzGitIt.github.io/sliR/reference/implied_allometry.md)),
`r_grad_app`, `r_grad_mass`, `b_anchor`, `beta_ref`, and its
decomposition `allometry_component` + `differential_component` (which
sum to `beta_ref` exactly, up to floating-point error).

## Details

The appendage-mass allometry is resolved exactly as in
[`implied_allometry()`](https://LezzGitIt.github.io/sliR/reference/implied_allometry.md)
(supply any two of `r_app_mass`, `b_ols`, `b_sma`, `b_avg`, or none for
the isometric default). The gradient enters through two further
correlations, `r_grad_app` and `r_grad_mass` –
[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)'s
own names for a gradient's correlation with `log(Append)` and
`log(Mass)` respectively.

Relative appendage size is defined against `b_anchor`: the scaling
exponent that counts as "no relative change". `beta_ref` is the
reference effect size a method anchored at `b_anchor` should report:

\$\$\beta\_{ref} = \rho\_{grad,app}\\ b\_{sma} - b\_{anchor}\\
\rho\_{grad,mass}\$\$

Writing \\\rho\_{grad,app} = \rho\_{grad,mass} + \Delta\rho\\, this
splits exactly into two components, returned separately because they
answer different questions about *why* `beta_ref` is non-zero:

\$\$\beta\_{ref} = \underbrace{\rho\_{grad,mass}(b\_{sma} -
b\_{anchor})}\_{\text{allometry\\component}} +
\underbrace{b\_{sma}\\\Delta\rho}\_{\text{differential\\component}}\$\$

`allometry_component` is non-zero whenever the realised allometry
departs from the anchor, even with no differential association between
the gradient and either trait (`r_grad_app = r_grad_mass`).
`differential_component` is non-zero only when the gradient associates
more strongly with one trait than the other. An anchor fixed a priori
(e.g. isometry, `b_anchor = 1/3`) is sensitive to both; anchoring at the
realised allometry itself (`b_anchor = b_sma`) zeroes
`allometry_component` identically, so `beta_ref` can only reflect
differential association. That is the distinction between anchoring
[`calc_sli()`](https://LezzGitIt.github.io/sliR/reference/calc_sli.md)
at a fixed exponent versus at each group's own estimated slope
(`control =`), and it changes what a significant `beta_ref` can be taken
as evidence of.

## See also

[`implied_allometry()`](https://LezzGitIt.github.io/sliR/reference/implied_allometry.md),
which this extends to a gradient;
[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)
to simulate data with this exact correlation structure.

## Examples

``` r
# Isometric anchor: both components can contribute.
implied_gradient_effect(b_sma = 0.5, r_app_mass = 0.3,
                         r_grad_app = 0.2, r_grad_mass = -0.1)
#> # A tibble: 1 × 12
#>   r_app_mass b_ols b_sma b_avg sd_log_append sd_log_mass r_grad_app r_grad_mass
#>        <dbl> <dbl> <dbl> <dbl>         <dbl>       <dbl>      <dbl>       <dbl>
#> 1        0.3  0.15   0.5 0.325         0.035        0.07        0.2        -0.1
#> # ℹ 4 more variables: b_anchor <dbl>, beta_ref <dbl>,
#> #   allometry_component <dbl>, differential_component <dbl>

# Anchored at the realised allometry itself: allometry_component vanishes,
# matching SLI-estimated's logic (calc_sli(control = ...)).
implied_gradient_effect(b_sma = 0.5, r_app_mass = 0.3,
                         r_grad_app = 0.2, r_grad_mass = -0.1, b_anchor = 0.5)
#> # A tibble: 1 × 12
#>   r_app_mass b_ols b_sma b_avg sd_log_append sd_log_mass r_grad_app r_grad_mass
#>        <dbl> <dbl> <dbl> <dbl>         <dbl>       <dbl>      <dbl>       <dbl>
#> 1        0.3  0.15   0.5 0.325         0.035        0.07        0.2        -0.1
#> # ℹ 4 more variables: b_anchor <dbl>, beta_ref <dbl>,
#> #   allometry_component <dbl>, differential_component <dbl>
```
