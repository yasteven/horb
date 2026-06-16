#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves plots

A0="${1:-1.5}"
MASS="${2:-1e11}"

TAG="dz2_xz_a${A0}_m${MASS}"
CSV="curves/${TAG}.csv"

echo "generating ${CSV}"
cargo run -q -p curve_fitter -- xz 3d_z2 "$A0" "$MASS" > "$CSV"

echo "plotting log density"
./scripts/plot_xz_density.py "$CSV" \
  --mode log \
  --contours \
  -o "plots/${TAG}_log.png" \
  --title "3d_z2 log density, a0=${A0} kpc, M=${MASS} M_sun"

echo "plotting linear normalized density"
./scripts/plot_xz_density.py "$CSV" \
  --mode linear \
  --contours \
  -o "plots/${TAG}_linear.png" \
  --title "3d_z2 normalized density, a0=${A0} kpc, M=${MASS} M_sun"

echo "plotting threshold/isodensity view"
./scripts/plot_xz_density.py "$CSV" \
  --mode threshold \
  --threshold 0.15 \
  --contours \
  -o "plots/${TAG}_threshold015.png" \
  --title "3d_z2 threshold density > 0.15 max, a0=${A0} kpc"

echo
echo "wrote:"
ls -lh \
  "plots/${TAG}_log.png" \
  "plots/${TAG}_linear.png" \
  "plots/${TAG}_threshold015.png"
