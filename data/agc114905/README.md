# AGC 114905 — 2024 Figure 7 scenario

Source: Mancera Piña et al. (2024), *Exploring the nature of dark matter with the
extreme galaxy AGC 114905*, [paper](https://arxiv.org/abs/2404.06537v2),
[DOI](https://doi.org/10.1051/0004-6361/202450230).

These are **digitized figure data**, not an author-supplied table. The five
circular-speed points, their error bars and the baryonic curves were extracted
from PDF vector paths in Figure 7 (PDF page index 8). The source PDF SHA-256,
coordinate calibration and original vector coordinates are retained in
`mp24_figure7_vectors.json`. `make_agc114905_data.py` reproduces the CSVs offline.

We hold the exact Figure 7 scenario fixed for every comparison: **D=80.27 Mpc,
i=35.86 degrees**. These are that paper's only-baryon fit parameters, not its
fiducial D=78.7 Mpc, i=31 degrees. No distance/inclination transformation or
additional geometry fit is applied. Radii are kpc, speeds and uncertainties km/s.
The points describe the circular speed after asymmetric-drift correction.

The red stars-plus-gas curve supplies `v2_baryons_kms2`. Its vector polyline is
interpolated as velocity and then squared; the original component polylines are
also exported for inspection. A 0.1 km/s digitization floor is added in quadrature
to the plotted errors. There is no published covariance supplied here, and
baryonic/distance/inclination uncertainties are not propagated. Residuals below
0.1 km/s should not be interpreted as accuracy beyond the extracted figure.

The 2024 data revise the earlier baryons-only interpretation: even the published
Figure 7 baryonic model falls short. In this extracted scenario, zero orbital
coefficients give **RMS 4.76 km/s, chi² 32.17 for five points**. The empty solution
is valid, but is not the best velocity fit. Four freely fitted weights plus
discrete subset selection can nearly interpolate these five points; that is an
in-sample flexibility result, not evidence for a uniquely recovered mass density.

```sh
./run.sh data agc114905
./run.sh --release do agc114905-fit          # reproduce recorded fits and graphs
./run.sh --release do agc114905-fit --refit  # repeat both complete subset searches
```

AGC 114905 is a separate dataset, not an extra member of the SPARC 175 sample.
