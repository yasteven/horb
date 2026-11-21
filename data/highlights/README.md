# Numerical records behind the four README figures

These small records are deliberately tracked along with the four PNGs. They
let readers inspect the exact parameters and plotted values without repeating
the full search. Intermediate force libraries and full candidate rankings
remain under ignored `reports/`.

For each of `mw` and `m31`:

- `*_fit.json`: selected superposition, original search settings and input
  hashes, convergence counts, masses, and the numeric fitted observation /
  baryon / signed-force snapshot.
- `*_weights.csv`: full-precision orbital weights used for the density export.
- `*_comparison.json`: fitted single-orbital and classical parameters, bounds,
  convergence flags, source links, scores and density-grid metadata.
- `*_scores.csv`: all six models on the same observations, including baryons.
- `*_curves.csv`: the six displayed total-speed curves, including values at
  every fitted observation radius.

Regenerate with `./run.sh --release do readme-highlights`; append `--refit` to
repeat both complete superposition searches. The [methods](../../docs/readme_highlights.md)
define “best”, the fit bounds, source assumptions and density color scales.
