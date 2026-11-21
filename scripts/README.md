# HORB plotting and analysis programs

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
- `measure_xz_lobe_extent.py`: morphology measurements.
- `fit_cuda_orbital_basis_wavelets.py`: existing positive density-basis fitter.
- `archive/make_mw_starter_curve.py`: historical starter data, available through
  `./run.sh python make_mw_starter_curve`.

Halo comparison plots and summaries include Einasto when present in the CSV.
Archived files without it still work and print a regeneration reminder.
The original fitter and baryonic prescription are unchanged.
Legacy shell workflows are available as `./run.sh do legacy-*` commands; the
original shell files remain recoverable through Git history.
