# sliR 0.2.0

Changes prompted by the first real use of the package, migrating the
SLI-allens-rule manuscript's simulation and analysis onto sliR.

* Added `implied_gradient_effect()`, extending `implied_allometry()` to a
  gradient: resolves the ground-truth reference effect (`beta_ref`) a
  correctly-specified estimator should recover, given an allometry plus the
  gradient's correlation with each trait (`r_grad_app`, `r_grad_mass`), and
  its two-pathways decomposition into an allometry-vs-anchor component and a
  differential-correlation component.

* Added `build_sli_slopes_hierarchical()`, and a `method =
  c("average", "hierarchical")` argument on `calc_sli()`. The hierarchical
  method resolves a group's slope from a reliability-gated cascade (own
  cell → single-variable marginal → pooled slope) rather than
  `build_sli_slopes_tbl()`'s unconditional per-variable averaging, falling
  back to a coarser grouping when a cell is too sparse or too weakly
  correlated to trust its own slope. Capped at two `control` variables.
  Unlike `build_sli_slopes_tbl()`, `NA`/`unknown_codes` rows are collapsed
  to an explicit `"Unk"` class and cascade like any other group rather than
  being dropped.

  Note: `calc_sli(control = NULL, method = "hierarchical")` and
  `calc_sli(control = character(0), method = "hierarchical")` are **not**
  equivalent — `control = NULL` is treated as "no control at all", falling
  back to the flat `b_sli` default, consistent with the `"average"` method.
  Constructions like `c(if (cond1) "Age", if (cond2) "Sex")` evaluate to
  `NULL`, not `character(0)`, when every condition is `FALSE` (`c(NULL,
  NULL)` is `NULL` in base R) — wrap in `as.character()` if you want that
  case to fall through to the pooled cascade instead.

* Added a "Getting started" vignette (`vignette("sliR")`) covering the two
  workflows not shown in the README: simulating across a parameter grid
  (including the gradient case), and computing the SLI per group on messy
  real data with unknown-coded and weakly-allometric groups.

* Added `sim_grid()`, a thin wrapper that maps [sim_allometric()] over a data
  frame of parameter combinations and row-binds the result, stamping a
  `.scenario` column identifying which row produced each simulated dataset.

* `sim_correlated()`'s `r` argument is renamed `r_app_mass`, matching
  `sim_allometric()`'s name for the same quantity. `r` still works but warns
  once per session; it will be removed in a future version.

* `inst/CITATION` now bakes the package version into the title, not just the
  `note` field: common CSL styles (APA, Elsevier-Harvard, Nature, Vancouver
  all checked) drop `note` for a Manual entry, which silently lost the
  version number that reproducibility depends on.

* Set up a pkgdown site (`_pkgdown.yml`, `.github/workflows/pkgdown.yaml`) so
  the reference and vignette render as a browsable site, deployed to
  `gh-pages` on push.

* `build_cov_mat()` gains an `sd_gradient` argument (default `1`). Previously the
  gradient always entered at unit standard deviation, so the returned matrix was
  a standardised structure rather than the covariance of the caller's actual
  system; a user had to rescale the gradient block by hand. The default is
  unchanged, so existing calls are unaffected.

* `?sim_allometric` now documents two behaviours that surprised a first user
  porting existing code: (1) the functions consume the random-number stream, and
  because the internal draw order is an implementation detail, exact
  reproducibility is tied to a pinned package version; (2) outlier trimming is
  applied *jointly* to appendage and mass, not sequentially, so the retained
  sample can differ slightly from code that trims one trait then the other.

# sliR 0.1.0

First release.

* `calc_sli()` computes the Standardized Length Index, adapting the scaled mass
  index of Peig & Green (2009) to standardise appendage length to a reference
  mass. Supports a scalar exponent or per-group SMA slopes via `control`.
* `build_sli_slopes_tbl()` and `build_group_cor_tbl()` fit and diagnose per-group
  allometric slopes.
* `implied_allometry()` resolves a full allometry from any sufficient subset of
  its parameters (the correlation and the OLS, SMA, and midpoint slopes).
* `sim_allometric()` simulates log-normal morphology with a target allometry,
  optionally along an environmental gradient drawn as a bounded, evenly-sampled
  design variable.
* `sim_correlated()` and `build_cov_mat()` provide the raw-scale bivariate draw
  and the underlying covariance matrix.
