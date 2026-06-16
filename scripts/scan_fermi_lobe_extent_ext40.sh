#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves reports plots

MASS="${1:-2.5e11}"
EXTENT="${2:-40}"
N="${3:-321}"

REPORT="reports/fermi_lobe_extent_scan_m${MASS}_ext${EXTENT}.csv"

echo "csv,threshold,z_top_kpc,z_bottom_kpc,height_total_kpc,x_at_top_kpc,x_at_bottom_kpc,r_at_top_kpc,r_at_bottom_kpc" > "$REPORT"

for a0 in 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.0; do
  CSV="curves/dz2_xz_a${a0}_m${MASS}_ext${EXTENT}.csv"

  echo "generating xz slice a0=${a0}, mass=${MASS}, extent=${EXTENT}"
  HORB_XZ_EXTENT_KPC="$EXTENT" HORB_XZ_N="$N" \
    cargo run -q -p curve_fitter -- xz 3d_z2 "$a0" "$MASS" > "$CSV"

  ./scripts/measure_xz_lobe_extent.py "$CSV" \
    --thresholds 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 \
    | tail -n +2 >> "$REPORT"

  for th in 0.15 0.20 0.25 0.30 0.40 0.50; do
    th_tag="$(printf "%.2f" "$th" | tr -d '.')"

    ./scripts/plot_xz_density.py "$CSV" \
      --mode threshold \
      --threshold "$th" \
      --contours \
      -o "plots/dz2_xz_a${a0}_m${MASS}_ext${EXTENT}_threshold${th_tag}.png" \
      --title "3d_z2 threshold >= ${th}, a0=${a0}, extent=${EXTENT} kpc"
  done
done

echo
echo "wrote $REPORT"
column -s, -t "$REPORT" || cat "$REPORT"
