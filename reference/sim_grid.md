# Simulate [`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md) across a grid of parameters

A thin wrapper for the common case of comparing many scenarios: each row
of `params` becomes one call to
[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md),
its columns matched to
[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)'s
arguments by name exactly as
[`purrr::pmap()`](https://purrr.tidyverse.org/reference/pmap.html)
would. Arguments that are constant across every scenario — most often
`gradient`, a single string naming one column — are passed via `...`
instead of appearing in `params`.

## Usage

``` r
sim_grid(params, ..., id_col = ".scenario")
```

## Arguments

- params:

  A data frame, one row per scenario, with columns named after
  [`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)'s
  arguments (e.g. `n`, `r_app_mass`, `r_grad_app`).

- ...:

  Arguments passed to every call of
  [`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md),
  constant across scenarios (e.g. `gradient = "Temperature"`).

- id_col:

  Name of the integer column identifying which row of `params` produced
  each simulated row, matching that row's position in `params`. Defaults
  to `".scenario"`.

## Value

A single tibble: the row-bound output of
[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md)
for every row of `params`, with `id_col` added. Join `params` back on by
row number (`dplyr::mutate(params, .scenario = dplyr::row_number())`) to
recover the parameter values behind each scenario.

## See also

[`sim_allometric()`](https://LezzGitIt.github.io/sliR/reference/sim_allometric.md),
which this maps over.

## Examples

``` r
grid <- expand.grid(n = 200, r_app_mass = c(0.1, 0.3, 0.5), b_sma = 1 / 3)
set.seed(1)
sims <- sim_grid(grid)
dplyr::count(sims, .scenario)
#> # A tibble: 3 × 2
#>   .scenario     n
#>       <int> <int>
#> 1         1   199
#> 2         2   200
#> 3         3   199

# Arguments constant across the grid, like the gradient's name and range,
# are passed once via `...` rather than repeated in every row of `params`.
grad_grid <- expand.grid(n = 200, r_grad_app = c(-0.3, 0, 0.3), r_grad_mass = -0.1)
set.seed(1)
grad_sims <- sim_grid(grad_grid, gradient = "Temperature", gradient_range = c(0, 24))
```
