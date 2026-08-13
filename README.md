
<!-- README.md is generated from README.Rmd. Please edit that file, then run rmarkdown::render("README.Rmd"). -->

# sliR

<!-- badges: start -->

[![R-CMD-check](https://github.com/LezzGitIt/sliR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LezzGitIt/sliR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

sliR computes the **Standardized Length Index (SLI)**, a size correction
for appendage lengths, and simulates allometric morphological data with
a known scaling exponent, correlation structure, and error.

The SLI adapts the scaled mass index of Peig and Green
([2009](#ref-peigNewPerspectivesEstimating2009)). Where their index
standardises *mass* to a reference *length*, the SLI standardises a
*length* to a reference *mass*:

$$SLI_i = A_i \left( \frac{M_0}{M_i} \right)^{b}$$

for individual $i$ with appendage length $A_i$ and body mass $M_i$, an
arbitrary reference mass $M_0$, and an allometric scaling exponent $b$.
The allometric scaling exponent is calculated using standardized major
axis (SMA) regression, as recommended by Warton et al.
([2012](#ref-wartonSmatr3Package2012)).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("LezzGitIt/sliR")
```

## Computing the SLI

Simulate a set of individuals, then size-correct their appendage lengths
under isometric scaling ($b = 1/3$, so linear dimensions scale with the
cube root of mass).

``` r
library(sliR)

set.seed(1)
birds <- sim_allometric(n = 500, b_sma = 1/3, r_app_mass = 0.3)

calc_sli(birds, b_sli = 0.33)
#> # A tibble: 498 × 5
#>    Append  Mass Append_log Mass_log   sli
#>     <dbl> <dbl>      <dbl>    <dbl> <dbl>
#>  1   178.  82.0       5.18     4.41  177.
#>  2   180.  78.5       5.19     4.36  181.
#>  3   174.  76.4       5.16     4.34  177.
#>  4   185.  77.0       5.22     4.34  188.
#>  5   184.  84.5       5.21     4.44  181.
#>  6   182.  90.3       5.20     4.50  175.
#>  7   178.  73.0       5.18     4.29  184.
#>  8   182.  77.5       5.20     4.35  184.
#>  9   185.  84.8       5.22     4.44  182.
#> 10   176.  75.7       5.17     4.33  180.
#> # ℹ 488 more rows
```

`calc_sli()` appends a column and returns the rows in their original
order. The function can be piped, creating several indices that sit side
by side for comparison.

``` r
library(dplyr)

est_b <- coef(smatr::sma(Append_log ~ Mass_log, data = birds))[["slope"]]

birds |>
  calc_sli(b_sli = 0.33,  rename_col = "sli_isometry") |>
  calc_sli(b_sli = est_b, rename_col = "sli_estimated") |>
  select(Append, Mass, sli_isometry, sli_estimated)
#> # A tibble: 498 × 4
#>    Append  Mass sli_isometry sli_estimated
#>     <dbl> <dbl>        <dbl>         <dbl>
#>  1   178.  82.0         177.          177.
#>  2   180.  78.5         181.          181.
#>  3   174.  76.4         177.          177.
#>  4   185.  77.0         188.          188.
#>  5   184.  84.5         181.          181.
#>  6   182.  90.3         175.          175.
#>  7   178.  73.0         184.          184.
#>  8   182.  77.5         184.          184.
#>  9   185.  84.8         182.          182.
#> 10   176.  75.7         180.          180.
#> # ℹ 488 more rows
```

$M_0$ is an arbitrary reference; changing it rescales every SLI by a
constant. The SLI maintains the units of the appendage $A$ and projects
it along the line implied by the scaling coefficient to the reference
mass $M_0$. The ranking of individuals will never change by changing
$M_0$. $M_0$ defaults to the mean of `Mass`.

If your columns are not named `Append` and `Mass`, you can name them
yourself — the arguments are unquoted:

``` r
lowercase <- rename(birds, wing = Append, mass = Mass)
calc_sli(lowercase, Append = wing, Mass = mass, b_sli = 0.33)
#> # A tibble: 498 × 5
#>     wing  mass Append_log Mass_log   sli
#>    <dbl> <dbl>      <dbl>    <dbl> <dbl>
#>  1  178.  82.0       5.18     4.41  177.
#>  2  180.  78.5       5.19     4.36  181.
#>  3  174.  76.4       5.16     4.34  177.
#>  4  185.  77.0       5.22     4.34  188.
#>  5  184.  84.5       5.21     4.44  181.
#>  6  182.  90.3       5.20     4.50  175.
#>  7  178.  73.0       5.18     4.29  184.
#>  8  182.  77.5       5.20     4.35  184.
#>  9  185.  84.8       5.22     4.44  182.
#> 10  176.  75.7       5.17     4.33  180.
#> # ℹ 488 more rows
```

## Group-specific allometric slopes

A single exponent for the whole sample assumes every group scales alike.
When age or sex classes scale differently, pass them to `control` and
each individual is corrected with their own group’s SMA slope.

``` r
birds$Sex <- rep(c("F", "M"), length.out = nrow(birds))
birds$Age <- rep(c("Juv", "Adult"), each = 2, length.out = nrow(birds))

build_sli_slopes_tbl(birds, control = c("Age", "Sex"))
#> # A tibble: 4 × 6
#>   Age   Sex       n b_sma_Age b_sma_Sex b_sli_avg
#>   <chr> <chr> <int>     <dbl>     <dbl>     <dbl>
#> 1 Adult F       124     0.304     0.317     0.310
#> 2 Adult M       124     0.304     0.340     0.322
#> 3 Juv   F       125     0.348     0.317     0.333
#> 4 Juv   M       125     0.348     0.340     0.344
```

Note the structure: with two control variables this fits *two* models —
mass interacted with age, and mass interacted with sex — and averages
their slopes per cell, rather than fitting a single age × sex
interaction. That keeps each slope estimated at the marginal rather than
the joint sample size.

Averaging hides disagreement, so when a cell’s per-variable slopes
differ by more than `slope_diff_warn` (0.15 by default) you get a
warning naming the cell with the largest discrepancy. Inspect the
`b_sma_*` columns before trusting `b_sli_avg`.

`calc_sli()` calls this for you. Individuals whose group is `NA` or
unknown-coded receive `sli = NA` and keep their row.

``` r
calc_sli(birds, control = c("Age", "Sex")) |>
  select(Append, Mass, Age, Sex, sli)
#> # A tibble: 498 × 5
#>    Append  Mass Age   Sex     sli
#>     <dbl> <dbl> <chr> <chr> <dbl>
#>  1   178.  82.0 Juv   F      177.
#>  2   180.  78.5 Juv   M      181.
#>  3   174.  76.4 Adult F      177.
#>  4   185.  77.0 Adult M      188.
#>  5   184.  84.5 Juv   F      181.
#>  6   182.  90.3 Juv   M      174.
#>  7   178.  73.0 Adult F      183.
#>  8   182.  77.5 Adult M      184.
#>  9   185.  84.8 Juv   F      182.
#> 10   176.  75.7 Juv   M      180.
#> # ℹ 488 more rows
```

Before trusting per-group slopes, check that each group *has* an
allometric relationship to fit. A group whose mass and appendage length
are uncorrelated yields a SMA slope that is noise.

``` r
build_group_cor_tbl(birds, control = c("Age", "Sex"))
#> # A tibble: 4 × 6
#>   Age   Sex       n b_ols     r    p_value
#>   <chr> <chr> <int> <dbl> <dbl>      <dbl>
#> 1 Adult F       124 0.626 0.409 0.00000236
#> 2 Adult M       124 0.201 0.144 0.111     
#> 3 Juv   F       125 0.421 0.319 0.000290  
#> 4 Juv   M       125 0.434 0.352 0.0000556
```

## Simulating allometric data

`sim_allometric()` draws from a bivariate log-normal specified by its
allometry rather than by raw standard deviations. That allometry has
only **two degrees of freedom**: any two of the correlation `r_app_mass`
and the three slopes `b_ols`, `b_sma`, `b_avg` determine the other two.
Supply two, or none to accept the isometric default
`b_sma = 1/3, r_app_mass = 0.3`.

The three slopes are not interchangeable. `b_ols` is the slope of
$E[\log A \mid \log M]$, and so the *generative* slope an `lm()`
recovers. `b_sma` is $\sigma_A / \sigma_M$, the slope an SMA fit returns
and the exponent `calc_sli()` assumes. Note that there is no clear
consensus regarding which line-fitting technique is best for allometric
scaling and related applications. Thus, `b_avg` chooses the midpoint
between `b_ols` and `b_sma`, which no estimator returns; it is a
convention for splitting the difference. `implied_allometry()` shows
what any parameter set implies:

``` r
implied_allometry(b_sma = 1/3, r_app_mass = 0.3)
#> # A tibble: 1 × 6
#>   r_app_mass b_ols b_sma b_avg sd_log_append sd_log_mass
#>        <dbl> <dbl> <dbl> <dbl>         <dbl>       <dbl>
#> 1        0.3   0.1 0.333 0.217        0.0233        0.07
```

``` r
set.seed(3)
d <- sim_allometric(n = 2000, b_sma = 1/3, r_app_mass = 0.3, trim_sd = NULL)

c(b_ols = coef(lm(Append_log ~ Mass_log, data = d))[[2]],
  b_sma = coef(smatr::sma(Append_log ~ Mass_log, data = d))[["slope"]])
#>     b_ols     b_sma 
#> 0.1000000 0.3333333
```

### Simulating along an environmental gradient

Pass `gradient` to add an environmental variable that both traits
respond to, and give its correlation with each log trait. **The gradient
is whatever you name it** — `"Temperature"`, `"Latitude"`, `"Year"`,
`"Rainfall"`, `"Elevation"`, `"Aridity"` — the string simply becomes the
column name. `gradient_range` then gives the units you want it expressed
in. Omit `gradient_range` and you get a z-scored gradient (mean 0, SD
1).

Because a gradient is a *sampling design* variable rather than a third
trait, it is drawn uniformly by default, so `gradient_range` is honoured
exactly and every part of the gradient is sampled equally.

The gradient and allometry arguments compose freely: the allometry
describes how appendage scales with mass, the gradient describes how
each responds to the environment.

``` r
set.seed(4)
gulls <- sim_allometric(
  n = 2000,
  # allometry: a shallower-than-isometric SMA slope, tightly correlated
  b_sma = 0.28, r_app_mass = 0.6, sd_log_mass = 0.09,
  # gradient: longer wings and heavier birds toward the poles
  gradient = "Latitude", gradient_range = c(-90, 90),
  r_grad_app = -0.3, r_grad_mass = 0.2,
  trim_sd = NULL
)

range(gulls$Latitude)
#> [1] -89.93252  89.93252
c(r_grad_app = cor(gulls$Latitude, gulls$Append_log),
  b_sma      = coef(smatr::sma(Append_log ~ Mass_log, data = gulls))[["slope"]])
#> r_grad_app      b_sma 
#>      -0.30       0.28
```

A normal gradient (`gradient_dist = "normal"`) reproduces the usual
three-variable multivariate normal, but is unbounded — it will happily
emit latitudes beyond $\pm 90$ — and leaves the extremes of the gradient
nearly unsampled.

Not every set of correlations is attainable. An appendage that lengthens
while mass falls along the gradient is incompatible with a strong
positive appendage–mass correlation, and `sim_allometric()` says so
rather than quietly approximating.

Note also that a correlation structure is **not** a causal structure. A
gradient correlated with both traits cannot distinguish a direct effect
on the appendage from an indirect one acting through body size.

### Measurement and biological error

Error comes in two flavours, each a fraction of the trait’s own standard
deviation. `meas_error` is observer error and hits both traits;
`transient_error_*` is genuine short-term biological fluctuation, which
afflicts mass far more than a skeletal appendage. Because that
fluctuation is independent of the appendage, it attenuates the observed
relationship:

``` r
set.seed(5)
clean <- sim_correlated(n = 2000, r_app_mass = 0.3)
set.seed(5)
noisy <- sim_correlated(n = 2000, r_app_mass = 0.3, transient_error_mass = 1)

c(clean = cor(clean$Append, clean$Mass),
  noisy = cor(noisy$Append, noisy$Mass))
#>    clean    noisy 
#> 0.300000 0.196199
```

`sim_correlated()` is the simpler raw-scale counterpart, with no log
transform and no allometric target — useful for building intuition about
how a fitted slope responds to correlation and to error on one trait but
not the other. `build_cov_mat()` exposes the covariance matrix that
`sim_allometric()` consumes.

## Function reference

| Function | Purpose |
|----|----|
| `calc_sli()` | Standardized Length Index per individual |
| `build_sli_slopes_tbl()` | Per-group SMA allometric slopes |
| `build_group_cor_tbl()` | Per-group diagnostic for whether allometry is meaningful |
| `sim_allometric()` | Log-normal traits with a target allometry, optionally along a gradient |
| `implied_allometry()` | What a set of allometric parameters implies for the remaining parameters |
| `sim_correlated()` | Raw-scale correlated appendage and mass |
| `sim_grid()` | Map `sim_allometric()` over a grid of parameter combinations |
| `build_cov_mat()` | Log-scale covariance matrix behind `sim_allometric()` |

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-peigNewPerspectivesEstimating2009" class="csl-entry">

Peig, Jordi, and Andy J. Green. 2009. “New Perspectives for Estimating
Body Condition from Mass/Length Data: The Scaled Mass Index as an
Alternative Method.” *Oikos* 118 (12): 1883–91.
<https://doi.org/10.1111/j.1600-0706.2009.17643.x>.

</div>

<div id="ref-wartonSmatr3Package2012" class="csl-entry">

Warton, David I., Remko A. Duursma, Daniel S. Falster, and Sara
Taskinen. 2012. “Smatr 3 – an R Package for Estimation and Inference
about Allometric Lines.” *Methods in Ecology and Evolution* 3 (2):
257–59. <https://doi.org/10.1111/j.2041-210X.2011.00153.x>.

</div>

</div>
