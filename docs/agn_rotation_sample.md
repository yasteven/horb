# HORB rotation curves, quasars and radio jets

Prepared 2026-11-25 from published observations and HORB's existing SPARC fits.

## The sample available now

There are **51 LeMMINGs–SPARC matches**, including 13 recovered through verified
NGC/UGC aliases. They provide 1,472 rotation measurements, activity measurements
or upper limits, and full-precision fitted orbital weights. Of these, 21 are
optical Seyfert/LINER nuclei. Twelve of those 21 have identified radio cores;
six further core detections occur in H II nuclei and are kept separate.
48 of the 51 pass the simple rotation-quality screen Q≤2 and inclination≥30°.
The screen is an analysis choice; all observations remain available.

The following seven are initial priorities based on resolved jet/lobe evidence.
Strength below is **published 1.5-GHz core luminosity**, in log10(W/Hz), and is
not an estimate of total mechanical energy output.

| Galaxy | SPARC name | Radio structure | log10 Lν | RC points | Outermost radius (kpc) |
|---|---|---|---:|---:|---:|
| NGC 1167 / B2 0258+35 | UGC02487 | Restarted compact jet, large outer lobes | 23.47 | 17 | 80.38 |
| NGC 2273 | UGC03546 | Twin jet plus extended emission | 20.46 | 30 | 29.23 |
| NGC 5005 | NGC5005 | Double-lobed structure, about 1.7 kpc | 20.11 | 18 | 11.47 |
| NGC 5985 | NGC5985 | Double-lobed structure | 19.58 | 33 | 34.72 |
| NGC 4051 | NGC4051 | Triple nuclear structure | 19.11 | 7 | 12.19 |
| NGC 7217 | UGC11914 | Twin lobes/triple source | 18.91 | 65 | 9.83 |
| NGC 2841 | NGC2841 | Twin jets/triple source | 18.66 | 50 | 63.64 |

The six small-scale morphology descriptions come from [LeMMINGs I,
Table A1](https://arxiv.org/abs/1802.02162); updated optical/radio measurements
come from [LeMMINGs III, Table 1](https://arxiv.org/abs/2109.06205).
For NGC 1167, see [Murthy et al. 2019](https://doi.org/10.1051/0004-6361/201935931)
and [Murthy et al. 2025](https://doi.org/10.1051/0004-6361/202453139).
Its jet power is estimated at roughly 10^44 erg/s using empirical radio-power
relations. The two published estimates, log10 Pjet=43.91 and 44.13, are
retained separately in `jet_power_estimates.csv`; they are not direct calorimetry.
Disturbed nuclear CO should not replace the extended disk's rotation curve.

The other five radio-detected optical AGN are NGC 2683, NGC 5033, NGC 2985
(UGC05253), NGC 3898 (UGC06787), and NGC 4395 (UGC07524). Their core/core-jet
classification is retained without upgrading every compact core to a resolved
jet. The full 51-object comparison group includes non-detections and nuclei
dominated optically by star formation. Lack of a LeMMINGs match for the other
124 SPARC objects is not evidence that they lack an AGN.

## Nearby quasar hosts

Yes: resolved host-galaxy rotation profiles exist for nearby quasars. The word
"closest" depends on whether a luminosity cut separates quasars from nearby
Seyfert nuclei. We preserve the source papers' sample classifications.

[Molina et al. (2021)](https://doi.org/10.3847/1538-4357/abd7f6)
observed six low-redshift PG hosts, obtaining rotation profiles for **five**:

| Target | z | Published luminosity distance (Mpc) | Rotation-data assessment |
|---|---:|---:|---|
| PG 0923+129 | 0.029 | 131.2 | Extended CO disk; outer flat part measured |
| PG 1011−040 | 0.058 | 267.9 | Profile measured, but kinematic disturbances; separate sensitivity sample |
| PG 1126−041 | 0.060 | 277.5 | Extended CO disk; outer flat part measured |
| PG 0050+124 / I Zw 1 | 0.061 | 282.3 | Extended CO profile plus a higher-resolution nuclear study |
| PG 2130+099 | 0.061 | 292.3 | Extended CO disk; outer flat part measured |

Their Figure 9 profiles reach at least about 8 kpc. The plotted velocities are
projected, so an importer must retain the documented inclination and geometry.
The public ALMA program is **2018.1.00006.S**, combined with ACA
**2017.1.00297.S**. [Fei et al. (2023)](https://doi.org/10.3847/1538-4357/acbb05)
provide a detailed I Zw 1 nuclear rotation/dispersion analysis inside about
2 kpc. Those nuclear and extended profiles overlap and must not be counted as
independent galaxies. PG 1244+026 is explicitly excluded from confident
kinematic analysis in the six-object study because it is too compact.

**Mrk 231**, z=0.04217, has a resolved rotating molecular disk, but its warped
nuclear disk, merger history and strong outflows require a more involved
analysis. [Feruglio et al. (2015)](https://doi.org/10.1051/0004-6361/201526020)
model the disk and outflow separately. It is useful for a central-dynamics
test, not yet a ready extended dark-halo curve in HORB.

The catalog stores each PG host's optical 5100-Å luminosity as a radiative
strength indicator. It does not infer a jet or jet power from that luminosity.
Numerical rotation points/errors have **not** been imported for these external
quasars. Obtain the published numerical profiles or measure them from the
archived cubes before fitting; the table identifies the exact next inputs.

## Powerful jets outside SPARC

**J2345−0449** is the first external priority: a spiral with a roughly 1.6-Mpc
double-double radio source, published Hα rotation out to about 25 kpc, and
independent ALMA CO-ring kinematics. Its total 1.4-GHz luminosity is about
2.5×10^24 W/Hz. Published empirical jet-power estimates are about
(1.7–2)×10^44 erg/s. Sources: [Bagchi et al. (2014)](https://arxiv.org/abs/1404.6889)
and [Nesvadba et al. (2021)](https://arxiv.org/abs/2103.12816).

[Ruffa et al. (2019), I](https://arxiv.org/abs/1901.07513) and
[II](https://arxiv.org/abs/1908.09229) provide six further radio galaxies with
resolved CO disks and kinematic modeling. Total radio luminosities in the
following table are comparable to one another; they have a different aperture
and frequency from the SPARC nuclear-core measurements above.

| Target | log10 total Lν at 1.4 GHz (W/Hz) | CO extent/role |
|---|---:|---|
| NGC 612 | 25.1 | About 9.6-kpc CO major-axis FWHM; also a much larger HI disk |
| IC 4296 | 24.2 | Nuclear disk, about 0.2-kpc major-axis FWHM |
| IC 1531 | 23.9 | Nuclear disk, about 0.25-kpc major-axis FWHM |
| NGC 7075 | 23.8 | Very compact disk; limited spatial resolution |
| NGC 3557 | 23.3 | Nuclear disk, about 0.3-kpc major-axis FWHM; SMBH contributes |
| NGC 3100 | 23.0 | About 1.6-kpc major-axis FWHM; warp and jet interaction |

The public ALMA program is **2015.1.01572.S**. These studies fit velocity cubes
and report circular-velocity models. Sampling their fitted arctan curves as if
they were independent observations would not be a valid new rotation test.
Their cube data and appropriate stellar/SMBH mass models must be used.

## Weight quantities prepared for the test

For signed component masses m_i, the export defines

`p_i = |m_i| / sum_j |m_j|`.

It records the largest p_i; effective component count `1/sum(p_i^2)`;
entropy `-sum(p_i log p_i)`; absolute negative-mass fraction; and cancellation
`1-|sum(m_i)|/sum(|m_i|)`. These are descriptors of the fitted density expansion,
not quantum probabilities. Raw signed weights, masses, scales and state labels
are also exported. A zero-mass expansion has blank concentration values.

The default dictionary and search remain exactly those of the SPARC reference
run: six m=0 orbital states, five scale multiples, up to four selected terms,
fixed stellar M/L, and the same signed-weight bound. These coefficients depend
on the chosen dictionary and baryonic assumptions. They are not a unique
reconstruction of the actual halo; global orbital masses can include material
well outside the measured rotation range. The export preserves radial coverage,
inclination, host luminosity, fit conditioning and coefficient-bound flags.

For the subsequent analysis, fix the sample and metrics before searching for
correlations. Compare detected radio-core strength, optical/[O III] strength
and X-ray strength separately, retaining censored measurements. Treat the seven
resolved jet/lobe cases as a morphology subset, not as the only observations.
Check dependence on host mass/luminosity, distance, radial coverage and inclination,
and repeat fits/bootstrap measurements to estimate how stable concentration is
under alternative near-equivalent fits. Present-day nuclear activity and older
outer lobes also trace different timescales. No correlation significance has
been claimed from this catalog preparation.

## Reproduction

```sh
./run.sh data agn
./run.sh check plots
# Optional fresh fits of all 51 matched galaxies:
./run.sh --release do sparc-fit --name agn51 --galaxies "$(cat reports/agn/galaxies.txt)" --jobs 2
./run.sh data agn --fits-dir reports/sparc/agn51 --out-dir reports/agn_refit
```

The default export already contains the existing reference fits; the fresh fit
commands are optional. New partial runs leave the other galaxies' fit metrics
blank rather than silently mixing their weights with the reference run.
See [data provenance and column conventions](../data/agn/README.md).


## Signed wavefunction lobes (2026-11-26)

The current comparisons also include a separate real-ψ gravitational source,
with intrinsic positive/negative lobes and absolute-mass normalization. See
[the complete model notes](signed_lobes.md).
`do signed-fit`, `do m31-signed-fit` and `do sparc-fit` now compare both definitions
by default; `--density-model squared` or `psi` selects one. SWELLS fits and projects
both families. `do agc114905-fit` handles the separate 2024 digitized galaxy dataset.
`do readme-highlights` regenerates all four galaxies' comparison figures.
`data agn --fits-dir reports/sparc/NAME --density-model psi --out-dir reports/agn_psi`
exports activity and ψ weights with absolute-strength concentration, not net-mass
concentration. `--max-components 0` explicitly selects zero orbital coefficients.

The signed-ψ activity export also works directly from the bundled 175-galaxy
fit records: `./run.sh data agn --density-model psi --out-dir reports/agn_psi`.
Use `--fits-dir` to select a different completed run.
