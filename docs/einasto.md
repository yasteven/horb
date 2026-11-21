# Einasto halo baseline

HORB implements a spherical Einasto halo alongside NFW, pseudo-isothermal and
Burkert. The convention is

\[
\rho(r)=\rho_{-2}\exp\left[-\frac{2}{\alpha}
\left(\left(\frac{r}{r_{-2}}\right)^\alpha-1\right)\right].
\]

`r_minus2` is the radius where the logarithmic density slope is −2;
`rho_minus2` is the density there. `alpha > 0` controls the curvature; another
common notation uses `n = 1/alpha`. Radii are in kpc, masses in solar masses,
densities in solar masses/kpc³, and circular speeds in km/s.

Writing `s=3/alpha`, `b=2/alpha`, and `t=b*(r/r_minus2)^alpha`, the mass is

\[
M(<r)=\frac{4\pi\rho_{-2}r_{-2}^{3}e^{b}}{\alpha b^{s}}
\gamma(s,t), \qquad v_c^2(r)=\frac{G M(<r)}{r},
\]

where `gamma` is the **lower** incomplete gamma function. The total mass at
infinite radius is finite. For the density, cumulative mass and total-mass
formulas, see Retana-Montenegro et al. (2012), *A&A* 540, A70, equations 4, 8,
12–13: [DOI](https://doi.org/10.1051/0004-6361/201118543),
[author-hosted paper](https://users.ugent.be/~mbaes/MyPapers/Retana-Montenegro%20et%20al.%202012a.pdf).

## Normalization and defaults

The constructor solves `rho_minus2` from `M(<r_ref) = m_ref`, just as the other
classical profiles do. This mass is **not** the Einasto mass at infinity.
HORB's orbital mass parameter continues to denote its total normalized mass;
the two conventions should be kept distinct when interpreting comparisons.

`alpha=0.18` is a reproducible starting value, not a universal constant or a
fitted result. In direct comparison commands, the Einasto radius defaults to
the supplied NFW scale radius. The MW workflows use 15 kpc by default.

`compare_dm_fixed_baselines` normalizes all four classical baselines to its
`baseline_mass` at `r_ref`, independently of the HORB mass. `compare_dm` uses
its single mass argument for the classical normalization as before.

## Commands

From the repository root:

```sh
# Standalone density, enclosed mass and circular velocity.
./run.sh run curve halo einasto 15 2.5e11 80 0.18

# Existing comparison arguments, optionally followed by Einasto radius and alpha.
./run.sh run curve compare_dm 3d_z2 1.0 2.5e11 5 15 8 80 15 0.18
./run.sh run curve compare_dm_fixed_baselines 3d_z2 1.0 2.5e11 2.5e11 5 15 8 80 15 0.18

# Generate comparison CSV, total-curve plot and scores.
./run.sh do mw-comparison
./run.sh do mw-comparison 1.0 2.5e11 2.5e11 20 0.25

# CUDA HORB versus the same four spherical classical baselines.
./run.sh do cuda-candidate 1.0 2.5e11 96 80 0.25 15.0 0.18
./run.sh graphs
./run.sh view
```

Old command forms remain valid. The optional standalone `halo` shape argument
is accepted only for Einasto. Nonpositive or nonfinite scale, shape, mass and
reference radius are rejected.

Default CPU comparison outputs share this stem:

`mw_lit_baryons_dm_compare_a1.0_m2.5e11_base2.5e11_einasto_rs15.0_alpha0.18`

The CSV is under `curves/`, the PNG under `plots/`, and the
`_score_inner.csv` / `_score_common.csv` tables under `reports/`.
Einasto parameters are included in workflow filenames so changing them does
not overwrite a differently parameterized comparison.

## CSV compatibility

The original first nine comparison columns are preserved. Five are appended:

| Column | Meaning |
| --- | --- |
| `v_einasto_kms` | Einasto dark-matter circular speed |
| `m_einasto_msun` | Enclosed Einasto mass |
| `einasto_r_minus2_kpc` | Radius parameter used for this curve |
| `einasto_alpha` | Shape parameter used for this curve |
| `einasto_rho_minus2_msun_per_kpc3` | Solved density normalization |

The parameters repeat on each row to make a CSV self-describing. The standalone
halo command retains its original four-column format. Plotters and summaries
accept archived comparison CSVs lacking Einasto and print a regeneration reminder.

## Comparison scores

All halo totals use the same existing literature baryon prescription, adding
velocity squares: `v_total² = v_halo² + v_baryons²`. This addition leaves the
baryonic model unchanged, including its existing spherical enclosed-mass
approximation for the exponential disk.

Each all-model report uses the same observations, errors and radius range.
The inner report covers 5–25 kpc; the common report requests 5–80 kpc, clipped
to the intersection of the sampled model domains. The CUDA comparison includes
both CUDA disk-plane HORB and spherical HORB. Observations outside the common
model domain are excluded rather than assigned the last sampled velocity.
The earlier `score_total_model_comparison_vs_mw` script uses this rule too;
its outer-range scores can therefore differ from archived reports that clamped
observations beyond 80 kpc. The separate Rust HORB-only scoring commands still
evaluate their requested radii directly.

Reports contain RMS residuals, chi-squared, and chi-squared **per point**.
The workflows evaluate fixed parameters; they do not optimize Einasto, assign
fitted degrees of freedom, or establish a best-fit ranking of halo families.
Signed orbital-density fitting is documented in [signed_density.md](signed_density.md).
SPARC/M31 ingestion remains later work.

## Numerical implementation and checks

`orbital_basis/src/halo_models.rs` evaluates the mass in log space. A convergent
lower-gamma series is combined analytically with its prefactor at small
arguments. Larger arguments use the regularized upper-gamma complement from
`statrs`, with default features disabled. This avoids overflowing Gamma or
underflowing the regularized lower gamma at small alpha/radius. As with the
other profiles, extreme parameter combinations must remain representable in
`f64`; failed normalization returns an error.

```sh
./run.sh check
./run.sh check plots
```

Checks cover high-precision reference values, an independent density integral,
the analytic `alpha=1` exponential case, central and infinite-radius limits,
the −2 slope, `dM/dr=4*pi*r²*rho`, mass scaling, CLI compatibility, invalid
parameters, and the comparison plotting/scoring programs. Python tests also
exercise archived CSVs and a narrower CUDA curve fixture. Executing the actual
CUDA solver requires the Jetson/GPU environment.
