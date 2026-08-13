# sliR — to do

Prioritised backlog, mostly surfaced by the first real use of the package
(migrating the SLI-allens-rule manuscript). Items done in the current dev
version are recorded in `NEWS.md`.

## High value

- [x] **Vignette for the two real workflows** (`vignettes/sliR.Rmd`, done
      2026-08-13). Covers simulating across a parameter grid (including the
      gradient case) and per-group SLI on messy real data (unknown-coded
      groups, gating on `build_group_cor_tbl()`).
- [ ] **README: reproducibility + citation section.** Recommend pinning a
      version (`@v0.1.0`) and citing a specific version/DOI, since outputs land
      in papers. Mint a Zenodo DOI from the GitHub release and reference it.

## Medium

- [ ] **Decide whether to re-add `calc_lambda()`.** It was dropped from sliR
      because nothing in the package used it, but the paper repos still do, and
      each keeps a private copy (`sqrt(calc_lambda(mass, wing))` approximates the
      SMA slope). Re-adding would let those repos delete the duplication.
- [ ] **`sim_correlated()` argument naming.** Its `r` is the same quantity as
      `sim_allometric()`'s `r_app_mass`; consider aligning for consistency
      (with deprecation of the old name).

## Lower / deferred (already noted in Project_notes / design)

- [ ] Grid-valued or integer gradients (e.g. `Year` as `1970:2026`) — needs a
      copula; currently the gradient is continuous.
- [ ] Multiple simultaneous gradients.
- [ ] The latent-size causal simulator (`gen_causal_data()` family) — high value,
      high maintenance; revisit on demonstrated demand.
- [ ] Stan hierarchical latent SEM — explicitly deferred.

## Housekeeping

- [ ] Set up pkgdown so the reference + vignettes render as a browsable site.
- [ ] Confirm the `inst/CITATION` version renders under common CSL styles
      (the Elsevier CSL dropped the `note = "R package version ..."` field).
