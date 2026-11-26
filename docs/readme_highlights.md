# README comparison figures

The MW/M31 figure pairs and J2141 rotation/lensing figures are generated together. Each rotation
plot shows the observations, their published errors, and seven models:

1. A fitted superposition of squared orbital densities, |ψ|².
2. An independently fitted signed-ψ lobe density, normalized by absolute mass.
3. Quantised inertia (QI), predicted from the fixed baryons without fitting.
4. Fitted NFW plus baryons.
5. Fitted Burkert plus baryons.
6. Fitted Einasto plus baryons.
7. Fixed baryons alone.

The single-orbital comparison remains available in the archived calculations,
but is omitted from the README figures. The three classical halos above share
the baryonic inputs with QI; see [the QI prescription](quantised_inertia.md).

## Reproduce

Install `requirements.txt`, select that Python with `HORB_PYTHON` if needed,
then run from the repository root:

```sh
./run.sh --release do readme-highlights
./run.sh check
./run.sh check plots
```

This regenerates the CPU force libraries, density slices, and MW/M31 figures
from the committed selected weights. It refits the three classical halos and
evaluates QI without fitting. It checks the numeric observation/baryon/force snapshot,
recomputes the fit scores, verifies that the weights match the report, and
independently checks the Rust density CSV against those weights before plotting.
Numeric comparisons tolerate roundoff between platforms; a changed dataset or
force calculation is rejected. No network access or CUDA is needed.

It also rebuilds J2141 from its recorded SWELLS fits: QI is projected through
the same observing model and scored on the six outer measurements. Lensing
selections are audited and redrawn without the single-orbital row. Recorded fits are reused. If the signed-ψ SWELLS record is absent, its rotation
fit is generated before predicting lensing; the squared solution is preserved.

```sh
./run.sh --release do readme-highlights --refit
```

`--refit` repeats the complete superposition searches and updates the recorded
winners. It can take several minutes. It does not commit or push; review changed
figures and records before publishing. `--dry-run` prints either command plan.

| Output | Location |
|---|---|
| MW/M31 README PNGs | `images/{mw,m31}_{rotation,density}_highlight.png` |
| J2141 README PNGs | `images/swells/j2141_{rotation,lensing}.png` |
| Vector PDFs | `plots/*_highlight.pdf` |
| Selected fit reports and exact weights | `data/highlights/{mw,m31}_{fit.json,weights.csv}` |
| Six MW/M31 curves and scores | `data/highlights/{mw,m31}_{curves.csv,scores.csv,comparison.json}` |
| J2141 observed predictions and scores | `data/highlights/j2141_{predictions.csv,comparison.json}` |
| Full search rankings on refit | `reports/{mw,m31}_showcase_candidates.csv` |
| Force libraries and density grids | `reports/*_showcase_basis.csv`, `curves/*_showcase_density_xz.csv` |

## Same observations, baryons and objective

The MW uses 29 Sofue (2020) observations in 5–80 kpc and the existing fixed
literature bulge + **sphericalized** exponential disk. The disk prescription
is retained from the earlier project comparison; it is not a thin-disk force.
Source: [Sofue (2020), Appendix A, Tables A1–A2](https://doi.org/10.3390/galaxies8020037).

M31 uses all 29 Corbelli et al. observations in 8.5–36.5 kpc with the fixed
Geehan Hernquist bulge, thin exponential disk and central point mass. It has
**no separately resolved gas component**. See [M31 provenance](../data/m31/README.md)
for the original sources, distance conventions and error-bar interpretation.
The fitted residual can absorb baryonic model error. No covariance matrices
are available; published error bars are used as diagonal weights.

Every fitted curve minimizes the same velocity-space statistic on its galaxy's
same observations:

\[
\chi^2 = \sum_j [(v_{\rm model}(R_j)-v_j)/\sigma_j]^2,
\qquad v_{\rm model}^2=v_b^2+v_{\rm DM}^2.
\]

The displayed RMS is the unweighted velocity residual RMS, not the optimized
statistic and not reduced χ². The score tables also record χ² and χ²/N. The
plots show only the fitted radial interval. Each README rotation figure uses
one full-height panel for all six curves, with per-model RMS values alongside
it. The separate signed-density figure follows below it in the README.
QI is a total-speed prediction from baryons, not an extra halo contribution;
its curve is scored on exactly the same observations without fitting.

## What “best superposition” means

The dictionary contains `1s,2s,2p_z,3s,3p_z,3d_z2` at each of
`a0 = 0.5,1,2,4,8 kpc`: 30 candidates. Every subset with zero through four
members is attempted, 31,931 in total. Weights are bounded to −20…+20, with
reference mass `1e11 M_sun`, and no weight normalization. For these showcase
figures the objective is **weighted χ², ridge=0**, rather than the default
fast workflow's sparsity/complexity score. The two workflows therefore need
not select the same winner.

The CPU calculation uses unsoftened axisymmetric multipoles, radial quadrature
order 128, 0.5-kpc sampling plus every observation radius, and interpolation
of signed force. Total circular force must remain nonnegative on all stored
nodes and three intermediate nodes per interval. The guard spans 0.5–100 kpc
for the MW and 0.5–60 kpc for M31. This is a sampled check, not a guarantee
outside that domain or at every unsampled point.

The recorded MW search has no failed subsets; M31 has two failed/infeasible
subsets, excluded from selection. All other candidates converged under the
existing constrained solver. Enumeration covers orbital choices on this scale
grid, not continuous scale values; nonlinear weight optimization is not a
proof of the global minimum. The selected weights are inside their bounds.

The MW winner has four positive terms. M31 has two positive and two negative
terms. Negative terms are used as signed density components, not quantum
probabilities. These flexible in-sample fits do not establish predictive
superiority or a unique physical dark-matter distribution.

## Baseline fitting and equations

NFW and Burkert fit mass and scale. Einasto fits mass, scale and α. Halo mass
means **enclosed mass at 80 kpc (MW) or 36.5 kpc (M31)**; it is not the total
mass of an infinite NFW or Burkert profile. Bounds are `1e8…1e13 M_sun`,
`0.05…500 kpc`, and `alpha=0.03…2`.

Fitting uses deterministic differential evolution (seed 42), followed by
bounded velocity-space least squares. Fit parameters, bounds and convergence
flags are recorded. These are refitted implementations of the published
profile families, not the original authors' published fitted parameter sets;
there are no concentration–mass or other cosmological priors.

The M31 Einasto fit reaches the lower shape bound `alpha=0.03`; its legend
identifies this. All other baseline parameters are interior to their bounds.
Thus this is a comparison over the stated bounded parameter domains, not an
unbounded optimum over every possible halo shape.

The classical force expressions match `orbital_basis/src/halo_models.rs`.
For NFW, enclosed mass is proportional to `ln(1+x) - x/(1+x)`; for Burkert,
`2 ln(1+x) + ln(1+x²) - 2 atan(x)`. Einasto uses the lower incomplete-gamma
ratio with shape `3/alpha` and argument `(2/alpha) x^alpha`. All are normalized
at the stated reference radius. Primary references:
[NFW](https://arxiv.org/abs/astro-ph/9611107),
[Burkert](https://arxiv.org/abs/astro-ph/9504041),
[Einasto mass integral](https://doi.org/10.1051/0004-6361/201118543).

The retained single-orbital implementation uses the exact `L=0,2,4` multipoles of `|R_32 Y_20|²`,
including exterior mass. Its radial density is proportional to
`x^4 exp(-2x/3)`; analytic incomplete-gamma integrals allow continuous scale
fitting without changing the nonspherical geometry. Tests compare it with
independent Rust multipole goldens and check the Newtonian far-field limit.
Tests also compare all three classical mass profiles with direct numerical
density integrals and recover synthetic fitted curves.

## Density rendering

Both 401×401 maps cover ±40 kpc with a ±12 kpc central zoom. They show the
**signed sum**, not its absolute value, in the `y=0` plane. Amber is positive,
cyan negative, and dark is zero. Symmetric-log normalization retains both
signs; each galaxy has its own labeled density scale (the zoom shares it).
Colors alone therefore should not be used to compare absolute densities
between galaxies. Baryons are not included in these dark-density maps.

The exact same weights produce the plotted force and density. The lower map
is a model implied by that choice of basis and coefficients, not an independent
observational measurement. Tests recompute all published scores from the
tracked plotted curves and check for mismatched weights or density grids.

Signed-ψ maps and records carry `_psi` in their filenames. AGC 114905 is also
regenerated by this command; its five digitized points and fixed Figure 7 geometry
are documented in [the data notes](../data/agc114905/README.md).
See [signed-lobe equations and normalization](signed_lobes.md).
