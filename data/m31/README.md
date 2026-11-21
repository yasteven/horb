# Andromeda (M31)

## Observations

`corbelli2010_table1.csv` contains all 29 rows and eight numeric columns of
Table 1 in Corbelli, Lorenzoni, Walterbos, Braun & Thilker (2010),
*A wide-field H I mosaic of Messier 31. II. The disk warp, rotation, and the
dark matter halo*, A&A 511, A89.

- [Paper and table](https://arxiv.org/html/0912.4133v1)
- [Journal DOI](https://doi.org/10.1051/0004-6361/200913297)
- [PDF](https://arxiv.org/pdf/0912.4133)

The CSV was extracted from the PDF text and every numeric cell was checked
against the paper's HTML table. This is tabulated data, not a traced graph or
the paper's fitted model. The radii are 8.5, 9.5, ..., 36.5 kpc. Speeds and
published error bars are retained verbatim numerically. In particular, the
large outer errors of 32.1, 32.0 and 38.2 km/s are not reduced or replaced.

The source adopts a distance of 785 kpc and a tilted-ring deprojection for
the warped HI disk. No second inclination correction is applied. The errors
follow the paper's Eq. 5: half the difference between the northern/southern
mean speeds, plus their mean dispersion. They are not standard errors of
the mean. We use diagonal weights because no covariance matrix is supplied.
The authors exclude the inner region because of non-circular motions.

M31 is absent from the [SPARC 175-galaxy catalog](https://astroweb.cwru.edu/SPARC/SPARC_Lelli2016c.mrt).
It is a separate dataset; the full SPARC importer is still future work.

## Fixed baryonic baseline

`geehan2006_model.json` records the rounded parameters from the
`(M/L_R)_d=3.3` column of Table 2 in
[Geehan et al. (2006), MNRAS 366, 996](https://arxiv.org/html/astro-ph/0501240v2).
These remain fixed while the orbital mixture is fitted.

| Component | Mass (solar masses) | Scale (kpc) | Force model |
| --- | ---: | ---: | --- |
| Bulge | 3.2e10 | 0.61 | Spherical Hernquist |
| Disk | 7.2e10 | 5.4 | Infinitesimally thin exponential disk |
| Central mass | 5.6e7 | — | Point mass |

The disk uses the Bessel-function force of Geehan Eq. 7,
`v_disk² = 2 G M_disk/Rd * y² * [I0(y)K0(y) - I1(y)K1(y)]`,
where `y=R/(2 Rd)`. This is the thin-disk force, not the old spherical
enclosed-mass approximation. The two source papers' distances (784 and
785 kpc) are retained without rescaling their published parameters.

This is a simplified literature baseline without a separately resolved gas
component. Gas and baryonic-model uncertainty can enter the fitted residual.
It is not a replication of Corbelli's fitted Sersic bulge, finite-thickness
disk and gas decomposition, and the resulting orbital density is not a
unique measurement of dark matter. The deprojected warped-disk curve is
compared to HORB's flat-plane axisymmetric force approximation.

## Run and reproduce

```sh
./run.sh data m31
./run.sh --release do m31-signed-fit cpu
./run.sh view plots/m31_signed_cpu_rotation.png
./run.sh view plots/m31_signed_cpu_density.png

./run.sh --release do m31-signed-fit cpu \
  --name m31_anneal --search anneal --seed 42 --budget 3000
./run.sh --release do m31-signed-fit cpu \
  --name m31_nonnegative --weights nonnegative
```

`do andromeda` aliases `do m31-signed-fit`. The menu's `a` entry runs the
CPU workflow. The `cuda` backend is available but needs GPU validation.

`data m31` works offline and creates `reports/m31_inputs/rotation_curve.csv`,
`baryons.csv` and `dataset.json`. The baryon table spans 0.25–100 kpc at
0.05 kpc spacing, includes each observation radius exactly, and stores
signed squared-speed force contributions. Interpolation is of force, not
speed, and extrapolation beyond its domain is rejected.

The fit defaults to 8.5–36.5 kpc, with force constraints on a 0.5–60 kpc
guard domain. The orbitals, scale grid and selection controls are the same
as the MW defaults. `--name` keeps separate experiments. `--rc-csv` and
`--baryons-csv` replace the inputs; automatic M31 provenance is omitted
when either is replaced. Supply matching `--dataset-meta` to retain a
custom provenance record. Metadata hashes must match the actual inputs.

Reports and figures use the `m31_signed_cpu` prefix. The JSON report embeds
source citations, data/model hashes, error conventions and model limitations.
`./run.sh dump --data` includes the bundled raw table and model parameters.
Generated inputs and outputs remain under ignored report/curve/plot folders.

## Example result

The default CPU run evaluated 504 subsets and selected:

| State | a0 (kpc) | Weight relative to 1e11 solar masses |
| --- | ---: | ---: |
| 1s | 8 | -1.58044185 |
| 2s | 2 | +2.23694082 |
| 2s | 8 | +12.0212697 |

The in-sample RMS is **3.432367 km/s**, with chi²/N **0.166422** over
29 measurements. These statistics depend on the fixed baryonic prescription
and the published error convention. They are not predictive validation or
a comparison to refitted classical halos. The data cover only 8.5–36.5 kpc;
the mixture's mass and behavior outside that interval are extrapolations.

Tests verify the tabulated measurements, error preservation, independent
high-precision thin-disk forces, interpolation of signed baryonic forces,
rejection of mismatched provenance, and M31-specific command routing.
