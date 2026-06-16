#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves plots reports

MASS="${1:-2.5e11}"

for a0 in 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.0; do
  CSV="curves/dz2_xz_a${a0}_m${MASS}.csv"

  echo
  echo "generating xz slice: a0=${a0}, mass=${MASS}"
  cargo run -q -p curve_fitter -- xz 3d_z2 "$a0" "$MASS" > "$CSV"

  for th in 0.05 0.10 0.15 0.20 0.25 0.30 0.40 0.50; do
    th_tag="$(printf "%.2f" "$th" | tr -d '.')"

    ./scripts/plot_xz_density.py "$CSV" \
      --mode threshold \
      --threshold "$th" \
      --contours \
      -o "plots/dz2_xz_a${a0}_m${MASS}_threshold${th_tag}.png" \
      --title "3d_z2 threshold rho/rho_max >= ${th}, a0=${a0}, M=${MASS}"
  done

  ./scripts/plot_xz_density.py "$CSV" \
    --mode log \
    --contours \
    -o "plots/dz2_xz_a${a0}_m${MASS}_log.png" \
    --title "3d_z2 log density, a0=${a0}, M=${MASS}"

  ./scripts/plot_xz_density.py "$CSV" \
    --mode linear \
    --contours \
    -o "plots/dz2_xz_a${a0}_m${MASS}_linear.png" \
    --title "3d_z2 normalized density, a0=${a0}, M=${MASS}"
done

echo
echo "regenerated plots:"
ls -lh plots/dz2_xz_a*_m${MASS}_threshold*.png | head
