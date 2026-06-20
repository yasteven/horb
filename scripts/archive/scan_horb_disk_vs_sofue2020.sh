#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

./scripts/make_mw_sofue2020_curve.py

MW="data/milky_way/sofue_2020_unified_rc.csv"
A0="${1:-1.5}"

REPORT="reports/score_horb_disk_scan_a${A0}_vs_sofue2020.csv"
echo "model_csv,r_min,r_max,n,rms_kms,chi2,chi2_per_point" > "$REPORT"

for HORB_MASS in 1.5e11 2e11 2.5e11 3e11; do
  for DISK_MASS in 6e10 8e10 1e11; do
    CSV="curves/dz2_total_a${A0}_m${HORB_MASS}_disk${DISK_MASS}_bulge1e10.csv"
    PNG="plots/dz2_total_a${A0}_m${HORB_MASS}_disk${DISK_MASS}_vs_sofue2020.png"

    echo
    echo "running a0=${A0}, M_HORB=${HORB_MASS}, M_DISK=${DISK_MASS}"

    cargo run -q -p curve_fitter -- total 3d_z2 "$A0" "$HORB_MASS" "$DISK_MASS" 3.0 1e10 0.7 > "$CSV"

    ./scripts/plot_total_vs_mw.py \
      "$CSV" \
      "$MW" \
      -o "$PNG" \
      --title "HORB a0=${A0}, M=${HORB_MASS}, disk=${DISK_MASS} vs Sofue 2020"

    ./scripts/score_total_vs_mw.py \
      "$CSV" \
      "$MW" \
      --r-min 5 \
      --r-max 25 \
      | tail -n 1 >> "$REPORT"
  done
done

echo
echo "inner scores, 5-25 kpc:"
column -s, -t "$REPORT" || cat "$REPORT"

echo
echo "best inner fits:"
sort -t, -k6,6n "$REPORT" | head -10
