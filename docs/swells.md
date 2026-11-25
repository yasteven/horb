# J2141 rotation fits and lensing predictions

This CPU workflow fits SDSS J2141−0001's gas rotation measurements with HORB,
NFW, Burkert and Einasto, freezes the mass parameters, and computes lensing
from those densities. Baryons alone and a single `3d_z2` orbital are included.
It compares with a published lensing summary; it does not yet reconstruct
the observed arc pixel by pixel.

## Run

```bash
./run.sh data swells
./run.sh do swells-fit
./run.sh view reports/swells/j2141/j2141_rotation.png
./run.sh view reports/swells/j2141/j2141_lensing.png

# Four HORB selections from bundled reference runs, plus all baselines.
./run.sh do swells-compare
./run.sh view reports/swells/j2141_comparison/j2141_lensing.png

# Reuse saved rotation parameters without another fit.
./run.sh do swells-fit --predict-only

# Nonnegative coefficient control, with the same data and baryons.
./run.sh do swells-fit --nonnegative --out-dir reports/swells/j2141_nonnegative

# Compare fresh matching runs instead of bundled reference runs.
./run.sh do swells-compare --signed-dir reports/swells/j2141 --positive-dir reports/swells/j2141_nonnegative

# Salpeter stellar-population sensitivity.
./run.sh do swells-fit --stellar-mass-offset-dex .24 --out-dir reports/swells/j2141_salpeter

# Optional exhaustive discrete search; much slower than the default beam.
./run.sh do swells-fit --exhaustive --budget 40000 --out-dir reports/swells/j2141_exhaustive

./run.sh check
./run.sh check plots
```

NumPy, SciPy and Matplotlib are existing dependencies. No CUDA is required.
The radial cumulative Simpson rule is implemented locally with NumPy, so
SWELLS does not require the `scipy.integrate.cumulative_simpson` function
introduced in SciPy 1.12. The supported older stack is Python 3.10, NumPy
1.21.5, SciPy 1.8 and Matplotlib 3.5.1. The same integration rule is used on
old and new SciPy; it is not replaced by lower-order trapezoidal integration.
Normal fitting commands write results; they do not commit each experiment.
The update installer checks and publishes explicitly enumerated source and
reference-result paths.

## Measurements and fixed inputs

Source: Dutton et al. (2011), [SWELLS II](https://arxiv.org/abs/1101.1622),
MNRAS 417, 1621. Table 4 supplies ten folded H-alpha/[N II] velocities and
errors. The **six outer points beyond 2 arcsec** are fitted, following the
source's central exclusion. All ten remain in the data and plots.

The default baryons use the independent Chabrier SPS masses in Table 2,
log10(M_bulge/M_sun)=10.26 and log10(M_disk/M_sun)=10.88, with the oblate
chameleon geometry in Table 3. The lensing-constrained masses are not inputs.
Cold gas is omitted because there is no resolved gas mass map in this source;
this is an assumption, not a measured zero. The inclination, stellar masses
and geometry are held fixed in each run. Uncertainties in these inputs have
not been marginalized.

The held-out Section 5.2 SIE summary is theta_E=0.89 (+0.05/−0.08) arcsec
and projected q=0.42 (+0.17/−0.12). Only the radius summary is compared
numerically. These are model-dependent arc summaries, not raw image data.

## Rotation implementation

The 30-term dictionary uses `1s,2s,2p_z,3s,3p_z,3d_z2` and scales
`0.5,1,2,4,8` kpc, with reference mass 1e11 M_sun. Up to four coefficients
in [-20,20] are fitted by velocity-space weighted chi-squared. Default beam
width is 8 and budget is 1200; the recorded run evaluates 636 subsets. This
is the best searched candidate, not a global optimum over arbitrary orbitals.

Analytic incomplete-gamma radial integrals retain all L=0,2,4 density
multipoles, matching the Rust force backend. Signed forces are added before
taking the square root. Component forces are never individually clipped.
Total inward equatorial force is guarded on 1e-4–10000 kpc grids, with
additional fine-grid constraints where cancellation reveals negative force
between coarse samples. This is not a test of full dynamical stability.

The observation operator applies inclination and flux-weights **velocity**
through a finite slit, spatial bin and Gaussian seeing PSF. It uses the
recorded dimensions and an exponential emitting disk. This approximation
does not reconstruct the actual dusty, ring-like emission distribution.
Every model has the same observing operator, errors, fitting mask and
baryonic components. Doubling the spatial/radial resolution audits its error.

NFW, Burkert and the single orbital fit mass and scale; Einasto also fits
alpha. Bounds are mass=1e8–1e13 M_sun, scale=0.05–500 kpc (single orbital:
0.05–50 kpc), alpha=0.03–2. Conventional halo mass means M(<13 kpc);
orbital mass means full-space mass. Differential evolution (seed 42) is
followed by least squares. RMS is unweighted and is not complexity-adjusted
model evidence; HORB also searches discrete orbital identities/scales.

## Lensing implementation

We assume the standard GR weak-field thin-lens relation, with the fitted
effective signed density sourcing both metric potentials equally. This is
a conditional mass-to-deflection prescription, not a derivation of a
relativistic HORB theory. See [Narayan & Bartelmann](https://arxiv.org/abs/astro-ph/9606001).

All angular-diameter distances are recomputed consistently with H0=70 and
Omega_m=0.3: 2.437736 kpc/arcsec and Sigma_crit=4.26769e9 M_sun/kpc². This
agrees with the source's Section 2 scale; the inconsistent distance values
printed in Section 4.4 are not mixed into the calculation.

Inclined m=0, l<=2 orbital densities project into exactly three angular
surface-density modes: cos(0 phi), cos(2 phi), cos(4 phi). Radial Poisson
integrals give the full two-dimensional deflection and Hessian, retaining
orbital shapes and signs without sphericalizing them or using an FFT box.
Classical halos are projected from their own spherical densities. Baryonic
chameleon deflections are analytic; density integrals set their mass
normalization consistently with the equatorial force.

For one closed tangential critical curve, theta_E,eff=sqrt(A/pi). Multiple
critical branches are flagged and **not assigned a single comparable radius**.
An outer-root envelope is stored only as a diagnostic. External convergence
and shear are zero in this pilot. The map exports reduced deflection,
convergence, shear and inverse magnification; no source reconstruction or
observed-image pixel chi-squared has been performed.

## Recorded results

Default fixed Chabrier baryons; same six rotation measurements:

| Model | LOS RMS (km/s) | Predicted effective Einstein radius (arcsec) |
|---|---:|---:|
| HORB signed superposition | 2.230 | 0.428 |
| HORB single orbital | 7.412 | 0.735 |
| NFW | 3.968 | 0.902 |
| Burkert | 4.271 | 0.838 |
| Einasto | 3.875 | 0.930 |
| Baryons alone | 93.880 | 0.630 |

The selected signed fit reaches its +20 coefficient bound. Its lensing radius
is below the published 0.81–0.94 arcsec interval; the three classical halos
lie within that interval. This preliminary comparison is not a likelihood
ratio or a model-selection result.

A separately selected **nonnegative-weight HORB control** gives RMS 3.009
km/s and radius **0.841 arcsec**. Other near-best signed subsets also produce
different lens predictions. Up to 20 alternatives with delta chi-squared<=1
are recorded, selected only by rotation. Their spread is not a confidence
interval; six rotation measurements do not determine a unique density.

### Separate rotation and lensing selections

`swells-compare` shows four HORB dots, with the same six rotation measurements
and fixed baryonic inputs. It defaults to the bundled reference runs:

| Coefficients | Selection | Rotation RMS (km/s) | Rotation χ² | Einstein radius (arcsec) |
|---|---|---:|---:|---:|
| Positive only | Best rotation | 3.009 | 5.216 | 0.841 |
| Positive only | Closest lensing | 2.990 | 5.502 | 0.898 |
| Signed (+/−) | Best rotation | 2.230 | 3.443 | 0.428 |
| Signed (+/−) | Closest lensing | 2.283 | 3.659 | 0.843 |

Best rotation means minimum **weighted χ²**, not minimum unweighted RMS.
Closest lensing means minimum absolute error from the published 0.89 arcsec
radius, among the **21 saved candidates** in each family: its rotation winner
and up to 20 alternatives with Δχ² ≤ 1. This is not a search of all possible
orbital weights, or of all 636/649 subsets evaluated during the original
signed/positive beam searches. Invalid or multiple-branch radii cannot win
the lensing selection. The observed radius is explicitly used to choose the
diamonds; these selected rows are not independent held-out lensing tests.
There is no continuous lensing refit of their rotation-fitted weights.

The density panel uses the **signed closest-lensing selection**, matching its
dot, critical curve and displayed rotation RMS. The command checks matching
inputs and baryons, reconstructs all four rotation scores, and recomputes
their lensing radii at 256 angular samples before drawing the figure.

The +0.24 dex stellar-mass sensitivity run gives HORB RMS about 1.81 km/s but
multiple critical branches, withholding a single-radius comparison. Thus
this pilot does not establish a unique HORB lensing prediction or a lensing
advantage over the conventional models.

## Outputs and checks

Each run writes `rotation_fit.json`, `rotation_predictions.csv`,
`lensing_predictions.json`, `summary.csv`, `surface_density_modes.csv`,
`critical_curves.csv`, `j2141_lens_map.npz`, and rotation/lensing PNG/PDFs.
`--fit-only`, `--predict-only` and `--no-plots` support staged runs.
Compact records are under `data/swells/reference/`; default figures are
under `images/swells/`.

`swells-compare` writes `selection.json` (including all selected weights and
selection provenance), `summary.csv`, `j2141_selected_lens_map.npz`, and the
updated lensing PNG/PDF to `reports/swells/j2141_comparison/`. The compact
comparison records are also bundled under `data/swells/reference/j2141_comparison/`.

Tests cover data transcription, distances, normalized forces and densities,
velocity-space projection, signed-weight recovery, analytic NFW deflection,
deflection derivatives and a known spherical critical radius. A separate
4320-point audit matches the Rust orbital backend to 3.1e-9 (km/s)^2.
Doubling lens projection, radial and angular resolution changes the default
Einstein radii by less than 1e-8 arcsec.

Five additional SWELLS objects are catalogued, awaiting original numerical
rotation arrays or documented digitization. They have not been fitted.
