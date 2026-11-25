# SWELLS spiral lenses

This extends HORB beyond the original SPARC 175. `catalog.csv` records six
spiral lenses with published gas rotation curves. **Only J2141 has numerical
rotation measurements imported and fitted here.** The other five entries
are a research inventory, awaiting numerical arrays.

- `j2141_rotation.csv`: all ten rows of Dutton et al. (2011),
  [SWELLS II](https://arxiv.org/abs/1101.1622), Table 4. These are observed
  line-of-sight velocities. The fitting mask excludes radii below 2 arcsec,
  following Section 6.4; the four inner points remain in the plots.
- `j2141_photometry.json`: independent SPS masses from Table 2, photometric
  chameleon geometry from Table 3, redshifts, cosmology and observing setup.
  The lensing-constrained masses from Table 5 are not used.
- `j2141_lensing.json`: the held-out, model-dependent SIE summary in Section
  5.2. This file is read only after rotation fitting.
- `reference/`: compact records of the signed fit, a nonnegative-weight
  control, and a +0.24 dex stellar-mass sensitivity run. These are conditional
  point fits, not marginalized parameter distributions.

These are attributed factual data transcriptions; no paper text or telescope
images are bundled. The other five curves appear in Figure 2 of
[SWELLS V](https://arxiv.org/abs/1206.4310). Its `J1331+3638` label is corrected
to `J1331+3628` in the [follow-up study](https://arxiv.org/abs/1609.05477).

```bash
./run.sh data swells
./run.sh do swells-fit
./run.sh view reports/swells/j2141/j2141_rotation.png
./run.sh view reports/swells/j2141/j2141_lensing.png
```

See [methods and recorded results](../../docs/swells.md).
