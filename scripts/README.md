# HORB plotting and analysis programs

## SWELLS spiral-lens comparison

`fit_swells.py` fits J2141's tabulated observed rotation velocities using a
slit/PSF observation model, then predicts lensing with the mass parameters
fixed. `horb_swells/` implements photometric baryons, signed force fitting and
non-spherical thin-lens projection. Run `./run.sh do swells-fit`; use
`./run.sh do swells-fit --help` for controls. See [methods](../docs/swells.md).
`compare_swells.py`, via `./run.sh do swells-compare`, draws four explicitly
labeled positive-only/signed selections. It ranks the saved candidates by
rotation χ² and by lensing-radius error separately, reconstructs their scores,
and uses the signed closest-lensing candidate for the density map.


Build/run/scan orchestration now lives in the Rust `tool/` crate. `run.sh` is
the only shell entry point. See [the tool guide](../tool/README.md) for every
workflow, the source dump command, dependencies and the shell migration table.

```sh
./run.sh do list
./run.sh python list
./run.sh plot xz-density --help
./run.sh python score_total_vs_mw --help
./run.sh python score_dm_comparison --help
./run.sh check plots
./run.sh --dry-run do cuda-basis
```

The Python programs can still be run directly. The tool invokes them with Python explicitly, so
executable file permissions and shebangs are not required.

- `plot_*.py`: density maps, rotation curves, halo comparisons and scan summaries.
- `score_*.py`, `rank_*.py`, `summarize_*.py`: diagnostics and rankings.
- `score_dm_comparison.py`: every halo with identical literature baryons and
  observations; optional CUDA curve; only scores the common sampled radius range.
- `horb_common/`: shared baryonic prescription, Einasto labels and comparison scoring.
- `tests/`: comparison, scoring and new/archived CSV plotting checks.
- `make_*contact_sheet*.py`: contact sheets from generated graphs.
- `make_mw_sofue2020_curve.py`, `convert_sofue2020_to_standard_rc.py`: MW input data.
- `make_m31_data.py`: offline Corbelli M31 curve, fixed Geehan thin-disk
  baryons and provenance with input hashes; invoked by `./run.sh data m31`.
- `measure_xz_lobe_extent.py`: morphology measurements.
- `fit_cuda_orbital_basis_wavelets.py`: existing positive density-basis fitter.
- `fit_signed_density.py`, `horb_signed/`: signed force-basis fitting, component
  search and constrained weight optimization; invoked by `./run.sh fit signed`.
- `plot_signed_fit.py`, `plot_signed_density.py`: fitted velocity, signed force
  contributions and signed density maps.
- `compare_signed_basis.py`: CPU/quadrature force-convergence checks.
- `archive/make_mw_starter_curve.py`: historical starter data, available through
  `./run.sh python make_mw_starter_curve`.

Halo comparison plots and summaries include Einasto when present in the CSV.
Archived files without it still work and print a regeneration reminder.
The original positive-only fitter remains available for historical workflows.
The new signed fitter preserves force signs. It accepts either the existing
component prescription or tabulated signed baryonic forces (`r_kpc`,
`v2_baryons_kms2`), including the M31 thin-disk model. Optional `--dataset-meta`
is verified against the actual input hashes and copied into the fit report.
Legacy shell workflows are available as `./run.sh do legacy-*` commands; the
original shell files remain recoverable through Git history.
## README figure generation

Use `./run.sh --release do readme-highlights` for all six featured figures.
`make_readme_highlights.py` validates and renders each matched pair;
`horb_highlights.py` supplies refitted classical forces and an independent
signed-density audit. `horb_qi.py` supplies the fixed QI prediction, while
`make_swells_highlight.py` rebuilds J2141 with its observing projection and
recorded lensing selections. Single-orbital curves/points are omitted. Source weights, plotted curves and score
tables are tracked in `data/highlights/`. See [methods and reproduction](../docs/readme_highlights.md).
## SPARC batch additions

Use `./run.sh data sparc` to normalize the bundled 175-galaxy tables, or
`./run.sh --release do sparc-fit --name sparc175 --jobs 2` to prepare, fit and
plot the full sample with verified resume. `make_sparc_data.py` and
`horb_sparc/` handle published inputs/provenance; `fit_sparc.py` orchestrates
the existing signed-density fitter and plotting programs. Details:
[SPARC workflow](../docs/sparc.md).

## AGN/jet catalog

`make_agn_catalog.py`, available as `./run.sh data agn`, joins verified SPARC
identifiers to the bundled nuclear radio and X-ray catalogs. It exports actual
rotation points and full-precision signed weights, with absolute-mass
concentration metrics. Radio luminosity, empirical mechanical jet power and
X-ray/optical radiative measurements retain their distinct units and meanings.
All source data and limits are attributed in `data/agn/`; the external quasar
and radio-galaxy target inventory records whether numerical rotation data have
actually been imported. See [the guide](../docs/agn_rotation_sample.md).
