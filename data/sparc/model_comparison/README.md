# SPARC model RMS and activity figures

These figures use the [paired 175-galaxy HORB snapshot](../density_comparison_20261126/README.md),
all 3,391 published SPARC measurements, and the [attributed activity catalog](../../agn/README.md).
Milky Way, M31, J2141 and AGC 114905 are separate datasets, not extra SPARC rows.

```sh
./run.sh python make_sparc_summary
./run.sh python make_sparc_summary --refit --jobs 4
```

The default command validates and reuses recorded classical fits. `--refit`
repeats all 700 classical halo fits; neither command reruns or modifies the
saved HORB search. Both write `images/sparc_model_rms.png` and
`images/sparc_jet_concentration.png`, with PDFs under `plots/`.

## Mean RMS

For each galaxy, RMS is `sqrt(mean((V_model - V_observed)**2))`, in km/s.
The plotted average is the arithmetic mean of these galaxy RMS values,
giving each galaxy equal weight. It is not pooled pointwise RMS, a median,
reduced chi-squared, or a test-set score. All fitted models minimize weighted
velocity chi-squared using the published errors; RMS is reported afterward.

Every model uses SPARC's fixed distance and inclination, gas forces, disk
M/L=0.5 and bulge M/L=0.7. Signed gas velocities retain their force signs.
All quality/inclination classes are included. The HORB reports' observation
and baryonic CSV hashes are verified against these exact inputs, for all 350
fits, before any aggregate is written.

| Model | Galaxies | Mean RMS (km/s) | Mean on common 174 (km/s) |
|---|---:|---:|---:|
| HORB squared density | 175 | 2.152 | 2.160 |
| HORB signed ψ | 175 | 2.691 | 2.695 |
| NFW | 175 | 7.487 | 7.507 |
| Burkert | 175 | 6.675 | 6.650 |
| Pseudo-isothermal | 175 | 6.231 | 6.214 |
| Einasto | 175 | 5.668 | 5.678 |
| VMS script profile | 175 | 10.515 | 10.528 |
| MOND, simple algebraic | 174 | 17.185 | 17.185 |
| QI | 174 | 19.428 | 19.428 |
| Baryons alone | 174 | 43.296 | 43.296 |

Full-precision values, including medians, are in `mean_rms.csv`.
**UGC01281** has baryonic v² = −33.0975 and −27.92125 (km/s)² at its first
two measured radii. Its QI, simple algebraic MOND and baryons-only full-curve
RMS values are undefined. No absolute value, zero-clipping, or deletion of
individual points is used to manufacture these three scores. `galaxy_scores.csv`
retains all ten model rows for all 175 galaxies, with an explicit status and
blank score for those three cases. White ticks in the figure show the means
on the common 174-galaxy intersection.

The two HORB variants use the saved beam search (width 6, budget 2,000,
up to four signed weights, |w|≤20, no ridge term). The squared dictionary has
30 columns and ψ has 20 identifiable columns; both search five scales and
states up to n=3. Each orbital/scale is a separate component. These are best
found candidates, not exhaustive searches or unique density determinations.

NFW, Burkert and pseudo-isothermal fit mass within the outermost measured
radius and scale radius. Einasto also fits alpha. Bounds are
1–10^14 solar masses, 0.001–10,000 kpc, and alpha=0.03–2. An explicit zero-halo
candidate is checked whenever baryons alone are valid. These broad mass
bounds include dwarfs; no Milky-Way minimum halo mass is imposed. Fitting uses
deterministic differential evolution (seed 42, population multiplier 12,
maximum 600 iterations, tolerance 10^-10), followed by velocity-space least
squares. All 700 global/local optimizations converged; 112 reach at least one
configured bound. Parameters, bound flags and zero-halo selections are saved
in `halo_fits.json`; cached scores are recalculated from the parameters before
plotting. Broad bounds are search conventions, not cosmological priors.

QI uses [McCulloch (2017), Eq. 8](https://doi.org/10.1007/s10509-017-3128-6)
with Theta=8.8×10^26 m and z=0. MOND is the fixed simple algebraic approximation
with a0=1.2×10^-10 m/s², no external field and no modified-Poisson field solver
([Famaey & McGaugh 2012, Eqs. 42–43](https://doi.org/10.12942/lrr-2012-10)).
Neither fits galaxy parameters. This is one specified MOND prescription,
not every MOND formulation.

The **VMS script profile** is the previously implemented
[downloadable script's](https://www.vms-institute.org/dark-matter/VMS_DARK_Matter_script.py)
`v_halo²=A*r/(r+s)`, with one fitted nonnegative amplitude and
`s=max(observed radius)/3` fixed. Source SHA-256:
`80d6b5301b17f89e5a56daf7458a8ce97e13b8768b95d1f0aa46a6271976e74f`.
It is compared using the same signed baryonic forces and original errors as
every other model. This label does not equate the script to other density
formulas in the VMS PDF, nor reproduce the source's entire baryon/error
pipeline. Its convex one-amplitude chi-squared fit is solved by a bracketed
derivative root or the allowed boundary. All parameters are in `vms_fits.json`.

These are in-sample results with unequal model flexibility. HORB's orbital
and scale search adds flexibility beyond the reported number of weights.
The mean RMS ranking alone is not model selection or evidence for a unique
physical interpretation.

## Activity and concentration

The upper panels include the **21 optical Seyfert/LINER nuclei** in the
LeMMINGs/SPARC overlap: 12 measured radio cores and nine upper limits.
The seven source-identified resolved jet/lobe hosts are highlighted and named.
The x coordinate is log10(1.5-GHz core Lν / W Hz^-1), rescaled to the SPARC
distance. It is a radio-activity proxy, **not mechanical jet power**. Only one
SPARC member has a mechanical-power estimate in the bundled catalog, so a
multi-galaxy jet-power correlation cannot be constructed from it.

The lower panels compare those seven resolved-jet hosts, the 14 remaining
optical AGN, and 29 HII-classified nuclei. The 124 galaxies without a catalog
match and one absorption-line (ALG) nucleus with uncertain AGN status are
excluded. HII spectra are a measured classification, not proof that an AGN or
jet is absent. This is therefore a descriptive comparison, not a secure
non-jet null sample. No luminosity limit or missing measurement becomes zero.

For each selected fit, concentration is

`C = max_i(|w_i * M_ref,i|) / sum_i(|w_i * M_ref,i|)`.

This uses **absolute component strength**, not integrated net mass or signed
probability. Both families have signed coefficients. The ψ normalization
uses the integral of |ψ|, so a zero-net-mass d orbital still contributes to C.
An all-zero expansion would have undefined concentration. The selection is
rotation-only; jet measurements never select the weights. Finite searches,
coefficient bounds and near-degenerate bases affect these numbers.

The resolved-jet medians are 84.7% (squared) and 82.6% (ψ), versus 74.6% and
66.5% in the HII group. Other optical AGN have medians 72.6% and 65.9%.
Among the **12 detected optical AGN cores only**, Spearman rho is −0.259
(p=0.417) for squared density and +0.252 (p=0.430) for ψ. Upper limits are
displayed but excluded from this ordinary rank calculation; these p-values
are descriptive, two-sided, and unadjusted. They do not establish a radio
luminosity correlation. No censored-data model or adjustment for stellar
mass, distance, morphology, rotation quality or selection effects is made.
No jet age, higher-n search, fit-uncertainty interval or independent replication
is included in this snapshot.

## Numerical records

- `galaxy_scores.csv`: 1,750 rows; scores, observation counts and undefined statuses.
- `mean_rms.csv`: equal-galaxy means, medians and common-sample means.
- `halo_fits.json`, `vms_fits.json`: fitted parameters and diagnostics.
- `activity_concentrations.csv`: both density families for all 51 measured
  activity matches, including the excluded ALG row, source-linked jet evidence,
  core limits, fit-bound flags, rotation scores and quality/inclination fields.
- `summary.json`: exact model/fit conventions, source hashes, coverage and
  activity statistics. The source catalog and selected weights are checked
  independently against each saved activity join before plotting.
