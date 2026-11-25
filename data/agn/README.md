# Activity measurements paired with SPARC rotation curves

This snapshot identifies **51 SPARC galaxies** in the 280-object LeMMINGs
survey. Their **1,472 measured rotation points** and existing HORB fits can
already be compared with nuclear activity measurements. There are 21 optical
Seyfert/LINER nuclei, 18 identified radio-core detections (12 in those optical
AGN), and 33 radio-core upper limits. Chandra supplies 31 nuclear X-ray
detections, 16 upper limits and four unobserved cases in the overlap.

```sh
./run.sh data agn
# To use a different completed SPARC fit run:
./run.sh data agn --fits-dir reports/sparc/sparc175 --out-dir reports/agn_latest
```

The default works offline, validates source hashes, and uses the full-precision
weights belonging to the committed SPARC reference run. It performs no new
curve fits or statistical correlation tests. See [the research and analysis
guide](../../docs/agn_rotation_sample.md) for target priorities and limitations.

| Bundled file | Contents |
|---|---|
| `lemmings_radio_2021.csv` | 280 nuclei: 1.5-GHz core/total luminosities, limits, optical class, radio morphology, [O III], Eddington-ratio estimates and published distances |
| `lemmings_radio_table1.tex` | Unchanged author Table 1; independently checked against the extracted radio luminosities and limit signs |
| `lemmings_xray_2022.csv` | 280 objects: 213 observed, with three X-ray bands, errors, limits, distances and Chandra observation identifiers |
| `simbad_aliases.csv` | Retrieved exact catalog identifiers for non-direct SPARC names; 13 additional matches |
| `jet_evidence.csv` | Seven priority SPARC objects with source-linked resolved jet/lobe evidence |
| `jet_power_estimates.csv` | Separate empirical mechanical-power estimates for NGC 1167 and J2345−0449 |
| `external_targets.csv` | Five PG quasar rotation-profile targets, Mrk 231, J2345−0449, six radio-galaxy CO disks, and one explicitly excluded compact PG source |
| `reference_fit_weights.csv` | 696 full-precision signed components across all 175 reference fits |
| `sources.json` | Provenance, URLs, extraction conventions and SHA-256 hashes |

Generated files under `reports/agn/`:

| File | Use |
|---|---|
| `sparc175_activity_weights.csv` | All 175 objects; the 124 without a catalog match remain unclassified |
| `matched51_activity_weights.csv` | The 51 measured/limited nuclei, ranked by detected core radio luminosity |
| `orbital_weights.csv` | Individual signed weights, masses, orbital scales and normalized absolute-mass fractions |
| `rotation_points.csv` | Published radii, velocities and errors, with signed baryonic force contributions |
| `galaxies.txt` | Exact comma-separated SPARC identifiers accepted by the fitter |
| `provenance.json` | Source/fit hashes, settings and counts for this export |

Radio luminosities in the original table are **νLν in erg/s**. Conversion to
Lν in W/Hz subtracts `7 + log10(1.5e9)` from their base-10 logarithms.
These monochromatic luminosities are **not mechanical jet powers**. Radio
core rankings are computed only for detections; upper limits never receive
ranks. Fixed power-decade labels describe radio luminosity, not an official
weak/strong AGN boundary. Published luminosities and distances are retained;
additional columns rescale luminosity by `(D_SPARC/D_published)^2` for distance
consistency with the rotation analysis.

Radio A means core/core-jet; B one-sided; C triple; D double-lobed; E complex.
An A or E classification alone does not establish a resolved jet. Even a
triple structure in an H II nucleus can have an uncertain origin; NGC 4217
is kept as `radio_origin_unresolved`. X-ray emission likewise can include
stellar or thermal contributions. An upper limit, missing coverage and an
optically star-forming classification are different conditions.

The physical observations are attributed third-party scientific data, not new
HORB measurements. Cite the original SPARC rotation references as well as:

- Baldi et al. (2018), LeMMINGs I, [10.1093/mnras/sty342](https://doi.org/10.1093/mnras/sty342).
- Baldi et al. (2021), LeMMINGs II, [10.1093/mnras/staa3519](https://doi.org/10.1093/mnras/staa3519).
- Baldi et al. (2021), LeMMINGs III, [10.1093/mnras/stab2613](https://doi.org/10.1093/mnras/stab2613), [author source, CC BY 4.0](https://arxiv.org/abs/2109.06205).
- Williams et al. (2022), LeMMINGs IV, [10.1093/mnras/stab3310](https://doi.org/10.1093/mnras/stab3310).
- [SIMBAD, CDS](https://simbad.cds.unistra.fr/), for alias resolution.

The 2026 [LeMMINGs VII 5-GHz release](https://arxiv.org/abs/2603.08562)
and its [public images](https://www.e-merlin.ac.uk/legacy-lemmings.html) are
identified for follow-up. Its full measurement table has not been imported in
this snapshot. The paper cautions that spectral indices need matched spatial
resolution and uv coverage; two observing frequencies do not by themselves
measure an intrinsic emission bandwidth.
