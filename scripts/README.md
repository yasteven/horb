# HORB scripts

## Active CUDA / Milky Way pipeline

- `make_mw_sofue2020_curve.py` — writes Sofue 2020 Milky Way rotation curve data.
- `convert_sofue2020_to_standard_rc.py` — converts Sofue data to standard RC schema.
- `scan_cuda_mw_3dz2.sh` — scans CUDA disk-plane `3d_z2` candidates.
- `scan_cuda_mw_3dz2_convergence.sh` — convergence scan over `64^3`, `96^3`, `128^3`.
- `run_cuda_mw_3dz2_candidate.sh` — plots one CUDA `3d_z2` candidate.
- `scan_cuda_mw_basis.sh` — scans multiple hydrogenic basis states.
- `plot_cuda_basis_scan_summary.py` — bar plot of best score per basis state.
- `plot_cuda_basis_bests_vs_mw.py` — overlays best basis-state candidates against Sofue 2020.
- `fit_cuda_orbital_basis_wavelets.py` — sparse positive density-basis fit.
- `scan_cuda_mw_wavefunction.sh` — scans true real wavefunction coefficients.
- `plot_cuda_wavefunction_scan.py` — plots wavefunction scan summaries.
- `rebuild_cuda_last_step_plots.sh` — rebuilds current CUDA/wavefunction plots.

## Active comparison / diagnostics

- `plot_cuda_total_vs_mw.py` — CUDA HORB + baryons vs classical halo models.
- `plot_cuda_vs_spherical_horb.py` — CUDA disk-plane force vs spherical approximation.
- `plot_dm_comparison.py` — DM-only HORB/NFW/pISO/Burkert comparison.
- `plot_total_model_comparison_vs_mw.py` — total model comparison against MW data.
- `plot_total_model_comparison_vs_standard_rc.py` — total model comparison against standard RC schema.
- `score_total_model_comparison_vs_mw.py` — scoring for total model comparisons.
- `summarize_dm_comparison.py` — summary of DM comparison CSVs.
- `summarize_total_curve.sh` — summary of total curve CSVs.

## Active Fermi / morphology pipeline

- `plot_xz_density.py` — XZ density slice plotter.
- `plot_orbital_views.sh` — log/linear/threshold views for one orbital candidate.
- `scan_morphology_scales.sh` — morphology scale scan.
- `regenerate_fermi_morphology_plots.sh` — rebuilds Fermi morphology plots.
- `regenerate_core_results.sh` — rebuilds core historical results.
- `scan_fermi_lobe_extent.sh` — Fermi lobe extent scan.
- `scan_fermi_lobe_extent_ext40.sh` — ext40 Fermi lobe extent scan.
- `rank_fermi_candidates.py` — ranks morphology candidates.
- `rank_fermi_lobe_extent.py` — ranks lobe extent candidates.
- `measure_xz_lobe_extent.py` — measures lobe extent from density CSV.
- `make_morphology_contact_sheet.py` — morphology contact sheet.
- `make_fermi_ext40_contact_sheet.py` — ext40 contact sheet.
- `make_fermi_log_linear_contact_sheets.py` — log/linear contact sheets.
- `make_fermi_threshold_contact_sheets.py` — threshold contact sheets.

## Archive

Old milestone, toy-baryon, starter-curve, and spherical-only scripts belong in `scripts/archive/`.
