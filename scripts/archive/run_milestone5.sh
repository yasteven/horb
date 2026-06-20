#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves plots reports

cargo fmt
cargo test -p orbital_basis

for a0 in 1.4 1.5 1.6 1.7; do
  csv="curves/compare_dm_horb_dz2_a${a0}_m1e11_vs_baselines.csv"
  png="plots/compare_dm_horb_dz2_a${a0}_m1e11_vs_baselines.png"
  report="reports/compare_dm_horb_dz2_a${a0}_m1e11_vs_baselines_summary.csv"

  echo
  echo "running HORB vs baselines for a0=${a0}"

  cargo run -q -p curve_fitter -- compare_dm 3d_z2 "$a0" 1e11 5.0 15.0 8.0 80.0 > "$csv"

  ./scripts/plot_dm_comparison.py "$csv" \
    -o "$png" \
    --title "HORB 3d_z2 a0=${a0} vs pseudo-isothermal, NFW, Burkert"

  ./scripts/summarize_dm_comparison.py "$csv" | tee "$report"
done

echo
echo "milestone 5 complete"
echo
ls -lh plots/compare_dm_horb_dz2_a*_m1e11_vs_baselines.png
