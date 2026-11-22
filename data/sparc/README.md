# SPARC source data

The two original machine-readable tables from
[the official SPARC site](https://astroweb.cwru.edu/SPARC/) are bundled unchanged
so every fit can run offline. Source URLs, retrieval date and SHA-256 checksums
are in [sources.json](sources.json). The loader verifies them before use.

**Citation:** Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016),
*SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate
Rotation Curves*, AJ 152, 157.
[DOI](https://doi.org/10.3847/0004-6256/152/6/157) ·
[Paper](https://arxiv.org/abs/1606.09251).
The catalog header retains the original rotation-curve reference keys and
their bibliography. Cite those observational sources when using individual galaxies.
These are attributed third-party scientific data, not newly authored HORB code.

| File | Contents |
| --- | --- |
| `SPARC_Lelli2016c.mrt` | 175 galaxy distances, inclinations, photometry, disk scales, quality flags and references |
| `MassModels_Lelli2016c.mrt` | 3,391 measured radii, velocities and errors, plus resolved gas/disk/bulge contributions |

All rows of the combined mass-model table were independently compared with the
175 individual files in the official `Rotmod_LTG.zip`; all numbers agree. The
combined table also agrees with the catalog distances. The ZIP is unnecessary
at runtime. Catalog quality counts are 99 high, 64 medium and 12 low; all remain
in the batch. Milky Way and M31 inputs remain separate from these 175 galaxies.

`Vgas` already includes helium. Disk and bulge velocities are normalized to
stellar mass-to-light ratio 1 at 3.6 microns. A negative component velocity
represents **outward radial force**, and its contribution is `V * abs(V)`;
it does not mean negative gas mass. Published errors exclude inclination
systematics. See [the fitting method](../../docs/sparc.md) for assumptions.

```sh
./run.sh data sparc
./run.sh --release do sparc-fit
./run.sh dump --data
```

Preparation exports normalized CSVs and per-galaxy provenance under
`reports/sparc_inputs/`. Batch fitting prepares its own matching inputs inside
each run directory. Source dumps with `--data` include both original tables.

## HORB reference run

[reference_summary.csv](reference_summary.csv) records all 175 completed default
CPU fits from 2026-11-21. [reference_run.json](reference_run.json) records the
settings, source/software hashes and validation counts. These two files are
HORB-generated results, separate from the unchanged observational tables.
Median in-sample RMS is 1.09 km/s; nine galaxies hit a coefficient bound, and
15 attempted subsets failed or were infeasible. All galaxies nevertheless have
an accepted fit. These diagnostics are retained for further refinement.
