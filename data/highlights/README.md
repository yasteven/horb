# Numerical records behind the README figures

These small records are deliberately tracked along with the PNGs. They
let readers inspect the exact parameters and plotted values without repeating
the full search. Intermediate force libraries and full candidate rankings
remain under ignored `reports/`.

For each of `mw` and `m31`:

- `*_fit.json`: selected superposition, original search settings and input
  hashes, convergence counts, masses, and the numeric fitted observation /
  baryon / signed-force snapshot.
- `*_weights.csv`: full-precision orbital weights used for the density export.
- `*_comparison.json`: fitted classical parameters, bounds, convergence flags,
  QI constants/redshift, source links, scores and density-grid metadata.
- `*_scores.csv`: all six models on the same observations, including baryons.
- `*_curves.csv`: the six displayed total-speed curves, including values at
  every fitted observation radius.

`j2141_predictions.csv` contains all ten observed LOS points, the six-point fit
mask and the six model predictions after the refined observing integration.
`j2141_comparison.json` records their scores, numerical audit, shared geometry,
baryons and QI redshift convention. All three rotation plots show HORB
superposition, QI, NFW, Burkert, Einasto and baryons alone.

Regenerate with `./run.sh --release do readme-highlights`; append `--refit` to
repeat both complete superposition searches. The [methods](../../docs/readme_highlights.md)
define “best”, the fit bounds, source assumptions and density color scales.
