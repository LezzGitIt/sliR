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

## Lower / deferred (already noted in Project_notes / design)

- [ ] Grid-valued or integer gradients (e.g. `Year` as `1970:2026`) — needs a
      copula; currently the gradient is continuous.
- [ ] Multiple simultaneous gradients.
- [ ] The latent-size causal simulator (`gen_causal_data()` family) — high value,
      high maintenance; revisit on demonstrated demand.
- [ ] Stan hierarchical latent SEM — explicitly deferred.

## Housekeeping

- [x] Set up pkgdown so the reference + vignettes render as a browsable site
      (done 2026-08-13; `_pkgdown.yml`, `.github/workflows/pkgdown.yaml`,
      `URL` field extended with the pages URL). Deploys to `gh-pages` on push
      to `main`/`master`; **GitHub Pages itself still needs to be pointed at
      the `gh-pages` branch in the repo's Settings → Pages** after the first
      push — not done here, as it changes shared repo settings.
- [x] Confirm the `inst/CITATION` version renders under common CSL styles
      (done 2026-08-13). Checked APA, Elsevier-Harvard, Nature, and Vancouver
      via `pandoc --citeproc`: all four drop the `note` field, so the version
      is now baked into the title as well, which survives every style tested.
