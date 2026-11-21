# HORB project operations

Run `./run.sh` from the repository root for the menu. The launcher also works
when invoked by absolute path from another directory. All relative data and
output paths are interpreted from the repository root.

The Rust `horb-tool` crate owns process execution, parameter sweeps, CSV report
aggregation, graph opening and source dumps. `run.sh` is the only shell script.
The existing Rust/CUDA models and Python plotting/numerical programs retain
their calculations. This change does not implement the planned signed-density
fitter, Einasto profile or additional galaxies.

## Everyday commands

```sh
./run.sh                         # menu
./run.sh help
./run.sh doctor                  # dependency diagnostics
./run.sh build                   # CPU crates, including the tool
./run.sh build all               # CUDA shared library, then all Rust crates
./run.sh check                   # tool format check and CPU workspace tests
./run.sh check all               # same, with CUDA tests (requires CUDA hardware)
./run.sh --release build all

./run.sh do make-curve
./run.sh do mw-comparison
./run.sh do cuda-candidate       # a0=1.0, mass=2.5e11, n=96, extent=80, eps=0.25
./run.sh do orbital-views 1.5 1e11
./run.sh graphs
./run.sh view                    # latest graph
./run.sh view 2                  # number from graphs
./run.sh view plots/dz2_xz_a1.5_m1e11_linear.png

./run.sh dump                    # current source -> unique dumps/cat_all_horb_*.txt
./run.sh dump --data             # also include input tables
```

`b`, `c`, `r`, and `d` alias `build`, `check`, `run`, and `dump`. `source` also
aliases `dump`. With redirected stdin, bare `./run.sh` prints help instead of
waiting for menu input.

## Source dump for a conversation

`./run.sh dump` reads current files from disk. It includes tracked source with
uncommitted edits and untracked, nonignored source. The header records the
branch, commit, working-tree status and file inventory. Rust, CUDA, Python, the
launcher, manifests, lockfiles, documentation and requirements are included.
Use `--data` to include input CSV/DAT/MRT/TXT tables under `data/`.

Generated output directories, build trees, previous dumps, images, binary
files, symlinks, `.env` files and private-key extensions are excluded. This is
a source selection rule, not a general secret scanner: check the file before
sharing if credentials were written inside source or documentation. For a
checkout without Git metadata the tool walks the source tree with the same
exclusions. Git-ignored source is omitted when Git metadata is present.

Each default filename contains a timestamp with nanoseconds, so repeated dumps
have different names. An explicit output is supported and must not exist:

```sh
./run.sh dump --data -o dumps/horb_for_review_01.txt
```

The command prints the full output path. Share that text file whenever the
current working tree is newer than the remote repository.

## Dependencies and CUDA

Install Rust/Cargo with `rustfmt`. The tool itself has no CUDA dependency.
The plotting programs use Python 3, NumPy, Matplotlib and Pillow:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
HORB_PYTHON="$PWD/.venv/bin/python" ./run.sh do orbital-views
```

The default Python executable is `python3`. Existing active virtualenvs work.
Plots use the noninteractive Matplotlib `Agg` backend, including over SSH.

CUDA commands call `nvcc` directly to build `libhorb_cuda.so` from
`horb_orbital_engine.cu`. They use `/usr/local/cuda/bin/nvcc` when installed there,
otherwise `nvcc` on PATH. `HORB_NVCC` can override the executable. The default
architecture remains `sm_87` for the Jetson Orin Nano; use `HORB_CUDA_ARCH` on
another GPU. The tool prepends the CUDA library directory to `LD_LIBRARY_PATH`
and preserves the inherited search paths. A nested workflow builds the shared
library once. A new command/menu action builds it again, so changed kernels are
picked up.

`view` uses `xdg-open` on Linux or `open` on macOS. `HORB_VIEWER` overrides the
viewer with a single executable path. Arguments are passed directly; these
environment values are not shell snippets. A desktop session/display is needed
to open a graphical viewer. Over headless SSH, use `graphs` and copy the listed
image to your desktop.

## Direct binaries, plotting and analysis

```sh
./run.sh run list
./run.sh run curve -- curve 3d_z2 1.5 1e11 > curves/example.csv
./run.sh run cuda write_horb_cuda_curve_basis -- 3d_z2 1.0 2.5e11 64 80 0.25
./run.sh python list
./run.sh plot xz-density curves/dz2_xz_a1.5_m1e11.csv --mode linear -o plots/example.png
./run.sh python score_total_vs_mw --help
./run.sh data mw
```

The direct commands retain the underlying program's arguments. `run curve`
without further arguments emits the default `3d_z2` curve. `run nbody` exposes
the existing placeholder integrator. Python program names accept underscores
or hyphens, with or without `.py`. The `plot` alias also accepts the name
without its `plot_` prefix. Numerical commands keep stdout available for CSV;
tool command traces go to stderr. When using shell redirection, create its
parent directory first, or use a `do` workflow which handles directories.

`--release` and `--dry-run` go **before** the top-level command:

```sh
./run.sh --dry-run --release do cuda-basis
HORB_A0_LIST=0.9,1.0 HORB_DM_MASS_LIST=2e11,2.5e11 ./run.sh do cuda-scan
```

Dry runs print the child commands and report actions without running them or
creating scientific outputs. The launcher may compile the tool itself.
Existing `HORB_*` scan overrides remain supported; they are passed to the
original programs. Child failures stop the workflow and return a nonzero exit
status. Redirected CSV output replaces its destination only after the child
process succeeds. Python programs and CUDA binaries that write their own files
retain their original output behavior.

## Migration from shell scripts

`./run.sh do list` lists every workflow and its positional arguments. For the
`do` workflows, the old script stem is also accepted, with or without `.sh`. Paths
such as `scripts/` are not part of the command name. Historical workflows are
explicitly labeled `legacy-*`; their scientific assumptions are preserved.

| Previous shell script stem | New command after `./run.sh` |
| --- | --- |
| `cat_scripts` | `dump` (now covers the whole project) |
| `cuda_kernels/cuda/build` | `build cuda` (also builds Rust CUDA binaries) |
| `scripts_make_curve` | `do make-curve` |
| `scan_dz2_density_axes` | `do density-axes` |
| `scan_dz2_fermi_scales` | `do fermi-scales` |
| `run_cuda_mw_3dz2_candidate` | `do cuda-candidate` |
| `scan_cuda_mw_3dz2` | `do cuda-scan` |
| `scan_cuda_mw_3dz2_convergence` | `do cuda-convergence` |
| `scan_cuda_mw_basis` | `do cuda-basis` |
| `scan_cuda_mw_wavefunction` | `do cuda-wavefunction` |
| `rebuild_cuda_last_step_plots` | `do cuda-rebuild-plots` |
| `run_mw_literature_baryons_model_comparison` | `do mw-comparison` |
| `summarize_total_curve` | `do summary CSV` |
| `plot_orbital_views` | `do orbital-views` |
| `plot_current_xz` | `do current-xz` |
| `plot_fermi_candidate_views` | `do fermi-views` |
| `scan_morphology_scales` | `do morphology-scan` |
| `regenerate_fermi_morphology_plots` | `do fermi-plots` |
| `scan_fermi_lobe_extent` | `do fermi-extent` |
| `scan_fermi_lobe_extent_ext40` | `do fermi-extent40` |
| `regenerate_core_results` | `do core-results` |
| `run_best_candidate` | `do legacy-best` |
| `run_best_candidate_total_model_comparison` | `do legacy-comparison` |
| `run_milestone3` | `do legacy-milestone3` |
| `run_milestone4` | `do legacy-milestone4` |
| `run_milestone5` | `do legacy-milestone5` |
| `run_mw_sofue2020_overlay` | `do legacy-sofue-overlay` |
| `run_mw_starter_overlay` | `do legacy-starter-overlay` |
| `scan_horb_disk_vs_sofue2020` | `do legacy-disk-scan` |
| `scan_horb_mass_vs_sofue2020` | `do legacy-mass-scan` |
| `scan_mw_literature_baryons` | `do legacy-literature-scan` |
| `scan_mw_multi_baryons` | `do legacy-multi-scan` |
| `scan_standard_rc_horb` | `do legacy-standard-scan` |
| `score_current_sofue_candidates` | `do legacy-score-current` |
| `score_sofue2020_ranges` | `do legacy-score-ranges` |

The old default curve scripts used a preexisting, unsupported argument form;
their tool replacements explicitly request `curve`. Historical workflows that
referenced moved archive scripts now invoke the appropriate programs directly.
CSV aggregation handles quoted fields, and rankings use column names instead
of positional shell sorting. A few historical graph filenames now use the full
candidate stem printed in the command trace.
