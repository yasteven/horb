#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves plots reports

cargo fmt
cargo test -p orbital_basis

for a0 in 1.4 1.5 1.6 1.7; do
  csv="curves/dz2_total_a${a0}_m1e11_disk6e10_bulge1e10.csv"
  png="plots/dz2_total_a${a0}_m1e11_disk6e10_bulge1e10.png"

  echo
  echo "running total curve for a0=${a0}"
  cargo run -q -p curve_fitter -- total 3d_z2 "$a0" 1e11 6e10 3.0 1e10 0.7 > "$csv"

  ./scripts/summarize_total_curve.sh "$csv"

  ./scripts/plot_total_curve.py "$csv" \
    -o "$png" \
    --title "HORB 3d_z2 + toy baryons, a0=${a0} kpc"
done

echo
echo "milestone 4 complete"
ls -lh plots/dz2_total_a*_m1e11_disk6e10_bulge1e10.png
