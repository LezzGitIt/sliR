# Per-cell SMA allometric slopes, with reliability-gated fallback

An alternative to
[`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md)
for when some `control` cells are too sparse or too weakly correlated to
fit their own SMA slope. Rather than averaging every `control`
variable's marginal slope unconditionally, this resolves one slope per
cell from a cascade, finest to coarsest, stopping at the first level
that clears both a minimum sample size and a reliable mass-appendage
correlation (`cor_min`, `cor_p_max`):

## Usage

``` r
build_sli_slopes_hierarchical(
  df,
  Append = Append,
  Mass = Mass,
  control = character(0),
  n_min_cell = 50,
  n_min_marginal = 100,
  cor_min = 0.3,
  cor_p_max = 0.05,
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

  Character vector of **0, 1, or 2** grouping columns. Each must be
  character or factor. Level 1 (the combined cell) is only attempted
  with exactly two; with `character(0)`, every row resolves to the
  pooled slope.

- n_min_cell:

  Minimum `n` for level 1 (the combined cell). Defaults to `50`.

- n_min_marginal:

  Minimum `n` for level 2 (a single-variable marginal). Defaults to
  `100`.

- cor_min, cor_p_max:

  A level passes only when its raw-scale Pearson correlation is at least
  `cor_min` (default `0.3`) with `p < cor_p_max` (default `0.05`).

- unknown_codes:

  Values in the `control` columns that mark an unknown group and so
  cannot be assigned a slope under `method = "average"`.

## Value

A tibble with one row per observed combination of `control`'s collapsed
classes (or one row overall, when `control` is empty): the control
columns, `pass_level` (`"cell"`, `"marginal"`, `"pooled"`, or `"none"`),
and `b_sli_avg` – the resolved slope, named to match
[`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md)'s
output so both feed
[`calc_sli()`](https://LezzGitIt.github.io/sliR/reference/calc_sli.md)`(control = ...)`
identically. `NA` wherever `pass_level` is `"none"`.

## Details

1.  **The cell itself** – the full cross of both `control` variables
    (only possible, and only attempted, with exactly two). Gated by
    `n_min_cell`.

2.  **Single-variable marginals** – each `control` variable's own slope,
    pooling across the other. Gated by `n_min_marginal`; averaged
    together when both pass, used unblended when only one does.

3.  **The pooled slope** – fit on every row of `df`, regardless of
    `control`. Reached when no finer level passes, and used for every
    cell when `control` is empty.

`df`'s pooled mass-appendage correlation gates the whole cascade: when
it is not itself reliable, every cell's slope is `NA`
(`pass_level = "none"`) rather than falling back to an unreliable pooled
fit. This keeps a caller's downstream SLI defined for either all of
`df`'s rows or none, matching the shared-N convention other
body-size-correction methods need for a fair comparison across them.

Unlike
[`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md)
and
[`build_group_cor_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_group_cor_tbl.md),
rows whose `control` values are `NA` or match `unknown_codes` are
**not** dropped: they are collapsed to an explicit `"Unk"` class per
variable, which cascades through the same three levels as every other
class – not a special case or an automatic fallback to pooled.

## See also

[`build_sli_slopes_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_sli_slopes_tbl.md),
the unconditional-averaging alternative this generalises;
[`build_group_cor_tbl()`](https://LezzGitIt.github.io/sliR/reference/build_group_cor_tbl.md)
for the same per-cell reliability diagnostic without the cascade.

## Examples

``` r
set.seed(1)
d <- sim_allometric(n = 500)
d$Age <- sample(c("Juvenile", "Adult", NA), nrow(d), replace = TRUE, prob = c(.1, .7, .2))
d$Sex <- sample(c("F", "M"), nrow(d), replace = TRUE)
build_sli_slopes_hierarchical(d, control = c("Age", "Sex"))
#> # A tibble: 6 × 4
#>   Age      Sex   pass_level b_sli_avg
#>   <chr>    <chr> <chr>          <dbl>
#> 1 Adult    F     cell           0.297
#> 2 Adult    M     cell           0.318
#> 3 Juvenile F     marginal       0.324
#> 4 Juvenile M     marginal       0.333
#> 5 NA       F     cell           0.383
#> 6 NA       M     marginal       0.345
```
