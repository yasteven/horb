#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves plots

cargo run -q -p curve_fitter -- xz 3d_z2 1.5 1e11 > curves/dz2_xz_a1.5_m1e11.csv
cargo run -q -p curve_fitter -- xz 3d_z2 1.7 1e11 > curves/dz2_xz_a1.7_m1e11.csv

./scripts/plot_xz_density.py curves/dz2_xz_a1.5_m1e11.csv \
  -o plots/dz2_xz_a1.5_m1e11.png \
  --title "3d_z2 density slice, a0 = 1.5 kpc, M_DM = 1e11 M_sun"

./scripts/plot_xz_density.py curves/dz2_xz_a1.7_m1e11.csv \
  -o plots/dz2_xz_a1.7_m1e11.png \
  --title "3d_z2 density slice, a0 = 1.7 kpc, M_DM = 1e11 M_sun"

echo
echo "wrote:"
ls -lh plots/dz2_xz_a1.5_m1e11.png plots/dz2_xz_a1.7_m1e11.png
