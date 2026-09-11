# Portfolio

Personal portfolio for actuarial and data science applications.
Static HTML, CSS and JavaScript. No build step, no framework, no dependencies.

## Local preview

```
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/serve.ps1 -Port 8777
```

Then open http://127.0.0.1:8777/

## Regenerating figures and data

Every chart in `assets/img/` and the aggregate file in `assets/data/` is produced
from the underlying model or dataset by one script. Nothing is hand-drawn.

```
Rscript scripts/make_figures.R
```

Requires R with `deSolve`, `ggplot2`, `dplyr` and `tidyr`.

Produces the SIR and SEIR curves, the cholera SEIRB model, the Kuramoto
synchronisation sweep, the Kaprekar triple counts, the two dengue figures, and
`assets/data/dengue.json` which drives the interactive chart on the page.

## Verification notes

The browser lab reimplements the SIR and SEIR systems in JavaScript with a
fourth-order Runge-Kutta integrator, independently of the R code. At the default
parameters the two agree:

| | Peak infected | Peak day | Attack rate |
| --- | --- | --- | --- |
| R, deSolve | 300.77 | 38.5 | 94% |
| JavaScript, RK4 | 301 | 39 | 94% |

The Kaprekar counts computed in R (1, 1, 1, 9, 5 for n = 1 to 5) reproduce the
values established independently in PARI/GP.

## Data

The Baguio dengue surveillance data was obtained from a health office and is
**not** committed here (see `.gitignore`). Only derived aggregates are published.
`scripts/make_figures.R` reads the raw files from a local path that will need
adjusting to run elsewhere.

Two traps in that dataset, both of which will silently produce wrong answers:

- The files are named as incidence per 1,000 but contain integer case counts.
- Several barangay names contain commas, so naive comma splitting undercounts.

## Structure

```
index.html                the whole site
assets/css/style.css      tokens, layout, light and dark themes
assets/js/app.js          theme, reveals, live solver, interactive chart
assets/img/               generated figures
assets/data/dengue.json   aggregates for the interactive chart
scripts/make_figures.R    regenerates every figure and the data file
scripts/serve.ps1         local static preview server
```
