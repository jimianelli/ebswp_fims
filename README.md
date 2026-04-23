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

## GitHub Pages
This repository is configured so `quarto render` writes the website to `docs/`.
For GitHub Pages, set the repository Pages source to:

- Branch: `main`
- Folder: `/docs`

After pushing the repo, the site home page will be `index.html`, and the report
will be available at `qmd/fims-implementation.html`.
