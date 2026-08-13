# sliR (development version)

Changes prompted by the first real use of the package, migrating the
SLI-allens-rule manuscript's simulation and analysis onto sliR.

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
