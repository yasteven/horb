# SPARC batch fitting

The SPARC workflow fits signed HORB density combinations to **all 175 galaxies**,
using all **3,391 published rotation-curve points**. It reuses the tested orbital
gravity and sparse fitter used for Milky Way and M31, while replacing their
analytic baryonic models with each SPARC galaxy's resolved gas, disk and bulge
force curves. The source tables ship with the repository; no download or login
is required to run a batch.

## Run

From the HORB root, with the existing Python dependencies installed:

```sh
python3 -m pip install -r requirements.txt
./run.sh check
./run.sh check plots
./run.sh --release do sparc-fit --name sparc175 --jobs 2
./run.sh view plots/sparc/sparc175/index.html
```

The tool builds the CPU `curve_fitter` once, then launches per-galaxy fits.
`--jobs 1` is the default to conserve Jetson memory. `--jobs 2` or more runs
independent galaxies concurrently; numerical-library threading is limited to
one thread per child. CUDA is not required or used by this batch workflow.
The HTML gallery opens locally and links to every rotation plot, signed-density
plot, fit report and weight table. On a headless host, copy both `plots/sparc/`
and `reports/sparc/` together to retain the relative links.

For a pilot, or to revisit particular galaxies:

```sh
./run.sh --release do sparc-fit --name pilot --galaxies NGC2403,DDO154,NGC7814 --jobs 2
./run.sh --release do sparc-fit --name pilot --jobs 2
./run.sh do sparc-fit --help
```

The second command expands the same run to all 175. Repeating a command resumes
it: completed galaxies are skipped only after all recorded artifact hashes,
input/source hashes, executable/software hashes and settings match. Missing,
altered or interrupted outputs are recomputed. One failed galaxy does not stop
the remaining fits; any requested failure produces a nonzero final exit code
and appears in the summary with its log path. Ctrl-C stops children; rerun to
resume. A per-run lock prevents simultaneous writers.

Use a new `--name` when changing fit or plot settings. `--force` repeats selected
galaxies with the same settings; it does not combine different settings in one
summary. The summary always lists all 175, distinguishing completed, pending,
failed and unselected galaxies. Normal scientific runs do not perform Git writes.

## Baryons and observations

All tabulated radii, observed velocities and errors are preserved. The baryonic
force proxy is SPARC's Eq. 2 before taking a square root:

\[
v^2_{\rm bar}(R)=V_{\rm gas}|V_{\rm gas}|
 +\Upsilon_{\rm disk}V_{\rm disk}|V_{\rm disk}|
 +\Upsilon_{\rm bulge}V_{\rm bulge}|V_{\rm bulge}|.
\]

Defaults are fixed **disk M/L = 0.5** and **bulge M/L = 0.7** in solar units at
3.6 microns; override with `--ml-disk` and `--ml-bulge`. The multiplier applies
to signed squared speed, not speed. The published gas contribution already
includes the 1.33 helium correction, so it is not multiplied again. Negative
component velocities, present in 48 galaxies, retain their outward-force sign.
The existing baryon interpolator interpolates signed squared speed and refuses
extrapolation beyond the first and last measured radii.

Published distances and inclinations are fixed; the already deprojected
rotation curves are not corrected again. Their uncertainties and quality flags
are recorded, but are not nuisance parameters in this first batch. Error bars
are used as diagonal velocity weights; no unprovided covariance matrix or
inclination systematic is invented. No quality or inclination cuts are applied.
The quality flag is visible in both summary and gallery.

## Orbital dictionary and search

For each galaxy the default dictionary has 30 columns: six states
`1s,2s,2p_z,3s,3p_z,3d_z2` at five values of `a0_star`:

```text
a0_star = [0.125, 0.25, 0.5, 1, 2] × published stellar disk scale length
```

`--scale-factors` changes those multiples; `--scales` supplies absolute kpc
values instead. `--states` selects distinct states from the six above.
The scale dictionary is fixed before fitting; individual scale values are
not continuously optimized. This definition accommodates dwarfs and giant
disks using the same dimensionless dictionary.

Every density column integrates to a common per-galaxy reference mass:

\[
M_{\rm ref}=R_{\max}\,[\max_R V_{\rm obs}(R)]^2/G,
\quad G=4.3009\times10^{-6}\ {\rm kpc}\,({\rm km/s})^2/M_\odot.
\]

This is a numerical normalization, **not a measured halo mass** or an additional
fitted parameter. `--ref-mass` can override it. Physical signed component masses
are `weight × Mref` and are saved alongside the weights. The default coefficient
bound is ±20, configurable by `--max-abs-weight`.

The default search is beam search, width 6, budget 2,000 subsets, seed 42,
up to four components, minimizing observed-velocity χ² with ridge penalty 0.
The effective component cap is `min(requested cap, N−1, dictionary size)` to
avoid allocating four coefficients to a four-point curve. Every selected
orbital identity and scale is also an adaptive choice: the cap does not make
such short-curve fits independent predictions.

CPU multipoles use radial quadrature order 128. The basis includes 129 uniform
radial grid nodes and all observed radii, confined to the observed domain.
The fitter adds intermediate positivity checks and fits
`sqrt(v2_baryons + sum(weight × v2_basis))` against the observed velocities.
The combined force must be positive on that guard grid. This does not require
positive dark density everywhere and gives no guarantee outside the domain.

For a more expensive search, different baryons, or a nonnegative control:

```sh
./run.sh --release do sparc-fit --name exhaustive --search exhaustive --budget 40000 --jobs 2
./run.sh --release do sparc-fit --name anneal --search anneal --budget 10000 --seed 43 --jobs 2
./run.sh --release do sparc-fit --name positive --weights nonnegative --jobs 2
./run.sh --release do sparc-fit --name disk06 --ml-disk 0.6 --jobs 2
./run.sh --release do sparc-fit --name fast_tables --plots none --jobs 2
```

The default 30-column, ≤4-term exhaustive search considers 31,931 subsets per
galaxy with enough observations. Optimizer failures and infeasible subsets are
counted. “Best” means the best accepted fit found by the requested search.
Beam/annealing searches are not exhaustive. χ²/N and RMS are **in-sample**;
χ²/N is not reduced χ², and these results alone do not establish a physical
interpretation or a unique density reconstruction.

## Output layout

| Path under the named run | Contents |
| --- | --- |
| `reports/sparc/NAME/run.json` | Global settings, source hashes, software fingerprints and selection |
| `reports/sparc/NAME/catalog.csv` | All 175 published catalog records plus radial coverage |
| `reports/sparc/NAME/summary.csv`, `summary.json` | Status, quality, fit scores, selected terms, convergence/bound diagnostics and masses |
| `reports/sparc/NAME/GALAXY/rotation_curve.csv`, `baryons.csv`, `dataset.json` | Exact observations, signed baryonic forces and provenance |
| `reports/sparc/NAME/GALAXY/basis.csv` | Signed orbital force library and numerical settings |
| `reports/sparc/NAME/GALAXY/fit.json`, `fit_*.csv` | Selected fit, every evaluated successful candidate, predictions and coefficients |
| `reports/sparc/NAME/GALAXY/density_xz.csv` | Density slice from the exact selected weights; extent ±last measured radius |
| `reports/sparc/NAME/GALAXY/status.json`, `run.log` | Completion hashes, timings and command/error log |
| `plots/sparc/NAME/GALAXY_rotation.png`, `GALAXY_density.png` | Matching rotation/force and signed-density figures |
| `plots/sparc/NAME/index.html`, `overview.png` | Searchable gallery and sample summary plot |

`--plots rotation` skips density generation; `--plots none` keeps fitting tables
and the HTML summary without PNG generation. `--density-bins` controls the
density slice resolution (default 201). Data preparation alone is available as
`./run.sh data sparc`; `./run.sh dump --data` includes the source tables for chat
snapshots. The existing Milky Way/Andromeda README highlights remain unchanged.

Source attribution and preserved original reference keys are documented in
[data/sparc/README.md](../data/sparc/README.md).

## Validated reference run (2026-11-21)

All 175 galaxies completed with the defaults above, producing 175 rotation
figures and 175 density figures. The median in-sample RMS was 1.09 km/s.
Nine galaxies reached a coefficient bound; 15 attempted subsets failed or
were infeasible. These flags remain visible in the records. Every selected
model had full column rank. These results use the current fixed baryons and
finite dictionary; they are an initial fitting baseline.

The run passed artifact-hash checks and independent recomputation of all
reported RMS/χ² values from its saved observations/predictions. All 175 cache
entries were reused on restart. Removing one generated plot caused only that
galaxy to rerun, restoring byte-identical output. All 50 CPU Rust tests and
32 Python tests passed; CUDA was not tested for this CPU workflow.

The compact [reference summary](../data/sparc/reference_summary.csv) and
[run settings](../data/sparc/reference_run.json) are committed. Large generated
force libraries, density grids and per-galaxy PNGs remain under the ignored
run directories and can be reproduced with the command above.
