#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves reports plots

MASS="${1:-2.5e11}"
REPORT="reports/fermi_lobe_extent_scan_m${MASS}.csv"

echo "csv,threshold,z_top_kpc,z_bottom_kpc,height_total_kpc,x_at_top_kpc,x_at_bottom_kpc,r_at_top_kpc,r_at_bottom_kpc" > "$REPORT"

for a0 in 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.0; do
  CSV="curves/dz2_xz_a${a0}_m${MASS}.csv"

  echo "generating xz slice a0=${a0}, mass=${MASS}"
  cargo run -q -p curve_fitter -- xz 3d_z2 "$a0" "$MASS" > "$CSV"

  ./scripts/measure_xz_lobe_extent.py "$CSV" \
    --thresholds 0.05 0.10 0.15 0.20 0.25 0.50 \
    | tail -n +2 >> "$REPORT"

  ./scripts/plot_xz_density.py "$CSV" \
    --mode threshold \
    --threshold 0.15 \
    --contours \
    -o "plots/dz2_xz_a${a0}_m${MASS}_threshold015.png" \
    --title "3d_z2 Fermi-lobe threshold, a0=${a0}, M=${MASS}"
done

echo
echo "wrote $REPORT"
column -s, -t "$REPORT" || cat "$REPORT"
