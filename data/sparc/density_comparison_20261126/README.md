# Paired HORB fits across SPARC 175

Both density families completed all 175 galaxies (3,391 rotation measurements),
with fixed SPARC disk/bulge M/L=0.5/0.7 and all quality flags retained.
The rotation-only beam search uses width 6, budget 2,000, at most four weights,
chi² selection, no ridge penalty, and coefficient bounds ±20. For small datasets
the number of terms is capped at N−1. The squared dictionary has 30 columns;
psi has 20 identifiable even-l, m=0 columns. Both use the same five scale
factors (0.125, 0.25, 0.5, 1, 2 times the stellar disk scale).

| In-sample RMS statistic | Squared density | Signed psi | Zero orbitals |
|---|---:|---:|---:|
| Median (km/s) | 1.093 | 1.406 | 41.020 |
| Mean (km/s) | 2.152 | 2.691 | 43.296 |

Signed psi has lower weighted chi² in 19/175 galaxies in this search. These
are best candidates found by a finite search, not exhaustive or predictive
cross-validation results. The two dictionaries differ structurally; the fits
do not establish unique physical orbital populations.

- `summary.csv`: both scores plus the explicit zero-orbital baryonic baseline.
- `orbital_weights.csv`: all 1,392 selected components, with density model,
  coefficient, absolute strength and net mass kept separate.
- `fit_reports.json`: both reports for every galaxy, input hashes, optimizer
  settings, constraints, failures and mass accounting.
- `run.json`: the completed run's settings, source hashes and software fingerprints.
- `agn_squared.csv`, `agn_psi.csv`: the 51 measured activity matches for each
  model separately. Radio upper limits and unknown activity remain identified.

Normalization and mode identifiability: [signed-lobe methods](../../../docs/signed_lobes.md).
The positive/negative lobe sums in the fit reports are component sums before
spatial overlap. Concentration uses absolute component strength, including
zero-net-mass multipoles; it is undefined for zero coefficients.

```sh
./run.sh --release do sparc-fit --name psi_models --density-model both --jobs 2
./run.sh view plots/sparc/psi_models/index.html
./run.sh data agn --fits-dir reports/sparc/psi_models --density-model psi --out-dir reports/agn_psi
```

The run generates per-galaxy comparison plots, one density map for each family,
all candidate rankings and a resumable gallery under `reports/` and `plots/`.
These bulky intermediate files are regenerated rather than added to Git.
