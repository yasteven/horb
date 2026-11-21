# Signed orbital-density fitting

## Definition

The model implemented here is a signed linear density expansion:

\[
\rho_{\rm DM}(\mathbf{x})=\sum_{i\in S} w_i\rho_i(\mathbf{x}),\qquad
\rho_i=M_{{\rm ref},i}|R_{n_i\ell_i}(r;a_i)|^2
|Y_{\ell_i m_i}^{\rm real}(\theta,\phi)|^2,\qquad w_i\in\mathbb{R}.
\]

Each orbital density integrates to its positive reference mass. Its multiplier
can be positive, zero or negative, and the weights do not have to sum to one.
Negative weights act as negative gravitational source density in this model.
The full-space signed mass of a component is `weight * ref_mass_msun`.
The report gives positive, negative and net mass separately.

This is the proposed signed-density model. Ordinary quantum superposition
instead forms `|sum(c_i psi_i)|²`, which includes interference terms and remains
nonnegative; changing a quantum amplitude's sign does not create negative
ordinary probability density. The coherent-wavefunction and older positive
density fit programs remain available separately. For the quantum convention,
see [MIT's introduction to quantum mechanics, position measurement](https://ocw.mit.edu/courses/22-02-introduction-to-applied-nuclear-physics-spring-2012/0456d26b0767e6aab3bace1e6f86d78b_MIT22_02S12_lec_ch2.pdf).

The chosen squared-orbital functions form a finite fitting dictionary; no
orthogonality or completeness claim is required. Degenerate or nearly
degenerate density/force combinations can occur. Matrix rank and condition
number are recorded, and fitted coefficients need not be uniquely identified
by a one-dimensional rotation curve.

## Linear gravity and the fitted quantity

Poisson gravity is linear in the source density. Thus potentials and forces
combine with the same signed weights. Define each disk-plane force basis by

\[
K_i(R)=R\,\partial_R\Phi_i(R,0),\qquad
q(R)=v_{\rm baryons}^2(R)+\sum_i w_i K_i(R),\qquad
v_{\rm model}(R)=\sqrt{q(R)}.
\]

`K_i` has units `(km/s)²`; it is retained even when negative. A positive
ring-like density can itself produce an outward force inside its ring, so a
speed-only table cannot reconstruct the force basis. The new library therefore
stores `v2_signed_kms2` directly. Old `v_cuda_kms` libraries are rejected by the
signed fitter and must be regenerated. No absolute values, weight squaring,
clipping of the dark target, or renormalization of the mixed density is used.

Only the total `q`, after baryons, must be nonnegative to give a real circular
speed. The optimizer checks this on all library radii, three intermediate
subdivisions per interval, and every fitted observation. It rejects invalid
solutions. Negative dark contributions and locally negative total dark density
remain allowed. Positivity between samples, outside this domain, orbital
stability and a dynamical origin for negative density are not established by
the fit.

See [Bovy, Dynamics and Astrophysics of Galaxies, §12.3](https://galaxiesbook.org/chapters/III-01.-Gravitation-in-Elliptical-Galaxies-and-Dark-Matter-Halos_3-Multipole-and-basis-function-expansions.html)
for linear density/potential expansions and the multipole formula.

## Force backends

The **CPU multipole** backend is the default. For the implemented real
harmonics (`l <= 2`), the azimuthal mean of `|Y_lm|²` contains only Legendre
orders `L=0,2,4`. Expanding it as `sum a_L P_L(cos(theta))`, the disk-plane force
is evaluated as

\[
K(R)=4\pi G M_{\rm ref}\sum_L \frac{a_L P_L(0)}{2L+1}
\left[\frac{L+1}{R^{L+1}}\int_0^R |R_{n\ell}(s)|^2s^{L+2}\,ds
-L R^L\int_R^\infty |R_{n\ell}(s)|^2s^{1-L}\,ds\right].
\]

The radial integrals use dimensionless Gauss-Legendre quadrature with a remote
exponential-tail cutoff; default radial order is 128. The angular expansion
is exact for the supported harmonics. This includes nonspherical disk-plane
forces; it is not the spherical `GM(<R)/R` shortcut. CPU forces are unsoftened.

Default states are `1s,2s,2p_z,3s,3p_z,3d_z2`, all axisymmetric. Named x/y p and
d orbitals are also supported, along with `n:l:m` for `1 <= n <= 10`, `l <= 2`.
For nonaxisymmetric states, rotation curves describe their **azimuthal mean**;
they do not describe exact circular orbits in the full nonaxisymmetric field.
The signed x-z density map shows the full selected real-orbital densities.
Pairs such as `2p_x` and `2p_y` have identical azimuthally averaged curves.

The **CUDA quadrature** backend evaluates the same mean densities by direct
summation of spherical quadrature cells. It retains signed force, accepts
signed cell masses, and returns an error on failed CUDA calls. It does not
renormalize the quadrature masses; their sum must agree with the analytic
reference mass within `1e-4`. Default resolution is `128 x 32 x 64` radial,
polar and azimuthal samples, with zero softening. Radial-force convergence can
require larger angular/radial orders, especially near dense shells. A nonzero
`--softening` deliberately changes the CUDA force model and is recorded.

## Search and objective

The program selects orbital identities and scales from the supplied library,
tries component counts from zero through `--max-components`, and optimizes
their real weights. The null model contains baryons only. Scale selection is
discrete; refine `--scales` to search a finer scale grid.

For each subset, the optimizer minimizes velocity-error-weighted chi-squared
plus a configurable small quadratic weight penalty:

\[
\chi^2=\sum_j\left[\frac{\sqrt{q(R_j)}-V_j}{\sigma_j}\right]^2,
\qquad J=\chi^2+\lambda\sum_i w_i^2.
\]

The squared-speed linear solution initializes the calculation; it is not the
final fitting objective. SLSQP uses analytic gradients, signed coefficient
bounds and linear total-force constraints (a tiny positive floor of 1e-8
`(km/s)²` keeps the square-root derivative finite during non-null fits).
Only converged feasible results
are retained. With nonnegative observed speeds this fixed-subset objective is
convex on positive total force; selecting the subset remains combinatorial.
[SciPy SLSQP documentation](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html).

| Control | Default | Meaning |
| --- | --- | --- |
| `--states` | six states listed above | Candidate orbital identities |
| `--scales` | `0.5,1,2,4,8` | Candidate a0 values in kpc |
| `--max-components` | `4` | Maximum number of selected terms |
| `--search` | `beam` | Beam, seeded annealing, or exhaustive subsets |
| `--beam-width` | `6` | Subsets retained at each size for expansion |
| `--budget` | `2000` | Maximum distinct subset evaluations |
| `--seed` | `42` | Annealing RNG seed |
| `--max-abs-weight` | `20` | Signed coefficient bounds |
| `--ridge` | `0.0001` | Weight-penalty coefficient lambda |
| `--criterion` | `bic` | Rank by J + k log(N); `chi2` uses J alone |
| `--weights` | `signed` | `nonnegative` provides a control fit |

The `bic` name denotes a BIC-style **heuristic complexity score**. Searching
the discrete dictionary adds model-selection freedom beyond the k weight
parameters; this score is not a calibrated evidence ratio or model probability.
Beam and annealing searches report the best explored subset, not a guarantee
of the global optimum. Exhaustive mode refuses a run whose combinations exceed
the budget. Scores are in-sample, use the supplied errors, and retain the
existing fixed baryonic prescription. Classical halo defaults from part 2
are not optimized competitors to this fit.

## Commands and outputs

Install/update the Python requirements in your existing environment or use:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
HORB_PYTHON="$PWD/.venv/bin/python" ./run.sh --release do signed-fit cpu
```

With dependencies on PATH:

```sh
./run.sh --release do signed-fit cpu
./run.sh view
./run.sh view plots/mw_signed_cpu_density.png
./run.sh --release do signed-fit cpu --name mw_anneal --search anneal --seed 42 --budget 3000
./run.sh --release do signed-fit cpu --name mw_positive_control --weights nonnegative
./run.sh do signed-fit --help
```

`--rc-csv` accepts another standard `R_kpc,Vobs_kms,eV_kms` table and
`--baryons-csv` its fixed component model, or a signed force table with columns
`r_kpc,v2_baryons_kms2`. Additional `v2_*_kms2` component columns must sum to
the total. Interpolation preserves force signs and rejects extrapolation.
`--dataset-meta` attaches source/model provenance only when its hashes match
both input files. [M31 has a bundled workflow](../data/m31/README.md) through
`do m31-signed-fit`; SPARC importing remains future work.
The complete workflow defaults to a basis domain of 0.5–100 kpc and a fit range
of 5–80 kpc. Target radii are included in the basis calculation. Standalone
fits interpolate **signed force** inside the library domain and reject
observations outside it. Use `--basis-r-min`, `--basis-r-max` and `--r-step`
to adjust the stored domain/sampling. Use a separate `--name` for each run
whose outputs you want to keep.

The default workflow produces:

| Path | Contents |
| --- | --- |
| `reports/mw_signed_cpu_basis.csv` | Versioned signed force library with state, scale, normalization and backend metadata |
| `reports/mw_signed_cpu.json` | Search controls, input hashes, errors, signed masses, rank/conditioning and chosen terms |
| `reports/mw_signed_cpu_weights.csv` | Full-precision selected signed coefficients and component masses |
| `reports/mw_signed_cpu_curve.csv` | Baryonic, dark, total and per-component signed force; total speed |
| `reports/mw_signed_cpu_observations.csv` | Observations, fit-selection flag and supported model predictions |
| `reports/mw_signed_cpu_candidates.csv` | Converged candidate subsets and their scores |
| `curves/mw_signed_cpu_density_xz.csv` | Signed x-z density slice |
| `plots/mw_signed_cpu_rotation.png` | Rotation curve and signed force contributions |
| `plots/mw_signed_cpu_density.png` | Diverging symmetric-log positive/negative density map |

`basis signed` generates a reusable library. `fit signed` searches an existing
library without recalculating gravity. `run curve signed_density` regenerates
a density slice from the saved weights. All support `--help`.

## Checks and example run

```sh
./run.sh check
./run.sh check plots
./run.sh check all                 # requires the CUDA/GPU environment
```

Tests cover signed density/force linearity, the analytic 1s mass profile,
high-precision 3d_z2 multipoles, outward force inside ring-like densities,
independent direct quadrature, synthetic recovery of known signed coefficients,
component-count selection, seeded search reproducibility, total-force
constraints, schema rejection, CSV reconstruction and plotting.

For a convergence comparison using the same states/scales/radius grid:

```sh
./run.sh --release basis signed cpu --out reports/reference_basis.csv \
  --rc-csv data/milky_way/sofue_2020_standard_rc.csv --n-radial 256
./run.sh python compare_signed_basis --reference reports/reference_basis.csv \
  --candidate reports/mw_signed_cpu_basis.csv --tolerance 0.000001
```

The same command can compare a CUDA library to the unsoftened CPU reference;
increase CUDA quadrature orders until its force errors are acceptable. A 1%
default per-basis tolerance is a diagnostic, not an error guarantee for a
large cancelling signed combination. Check the final weighted curve too.

The default CPU example for this update evaluated 500 subsets and selected:

| State | a0 (kpc) | Weight, relative to 1e11 solar masses |
| --- | ---: | ---: |
| `3s` | 2 | −4.24361572 |
| `3p_z` | 2 | +9.52168934 |
| `3d_z2` | 1 | +0.42200467 |

On the 29 Sofue observations in 5–80 kpc, RMS is 4.430804 km/s and chi²/N is
0.103321. Positive mass is 9.94369401e11 solar masses, negative mass is
−4.24361572e11, and net mass is 5.70007829e11. Increasing radial order from 128
to 256 changed the basis force norms by at most 1.2e-14 relatively and the
selected weighted dark force by about 1.2e-9 (km/s)². This verifies numerical
convergence of this CPU example; it does not establish an independent
prediction or the physical existence of negative density. The CUDA numerical
execution must be checked on the Jetson/GPU.
