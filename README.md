# ebswp_fims

EBS walleye pollock FIMS repository with a Quarto website report.

## Repository layout
- `qmd/`: Quarto report sources.
- `scripts/`: data build and model-fit scripts.
- `data/`: derived model inputs.
- `outputs/`: model-fit artifacts used by the report.
- `docs/`: rendered website for GitHub Pages.

## Local workflow
1. Build data:
   - `Rscript scripts/01_build_data.R`
2. Fit model:
   - `Rscript scripts/02_fit_model.R`
3. Render the site:
   - `quarto render`

The main report source is `qmd/fims-implementation.qmd`.

## SparseNUTS framework

This repository now includes a first-pass framework for Bayesian sampling with
`SparseNUTS` using the existing FIMS/TMB objective.

Files:
- `scripts/sparsenuts_framework.R`: shared helpers to build the FIMS input,
  fit the conditional mode, extract the TMB object, and run `SparseNUTS`.
- `scripts/03_fit_sparse_nuts.R`: command-line runner for the SparseNUTS path.

Example:

```sh
Rscript scripts/03_fit_sparse_nuts.R \
  --num-samples=250 \
  --num-warmup=250 \
  --chains=4 \
  --cores=1 \
  --metric=diag \
  --seed=123
```

This writes `outputs/sparsenuts_fit.rds` by default.

Notes:
- The current framework reuses the existing `fit_fims()` mode-finding step as
  the preconditioning stage for SparseNUTS.
- It operates on the TMB object created by the FIMS application rather than a
  separate RTMB rewrite.
- `cores=1` is the safest default for initial debugging. Increase parallelism
  only after confirming the model runs cleanly in serial.

## GitHub Pages
This repository is configured so `quarto render` writes the website to `docs/`.
For GitHub Pages, set the repository Pages source to:

- Branch: `main`
- Folder: `/docs`

After pushing the repo, the site home page will be `index.html`, and the report
will be available at `qmd/fims-implementation.html`.
