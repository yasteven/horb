#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

./scripts/make_mw_sofue2020_curve.py

MW="data/milky_way/sofue_2020_unified_rc.csv"
REPORT="reports/score_horb_total_vs_sofue2020.csv"

echo "model_csv,r_min,r_max,n,rms_kms,chi2,chi2_per_point" > "$REPORT"

for a0 in 1.4 1.5 1.6 1.7; do
  for horb_mass in 4e10 5e10 6e10 7e10; do
    csv="curves/dz2_total_a${a0}_m${horb_mass}_disk6e10_bulge1e10.csv"
    png="plots/dz2_total_a${a0}_m${horb_mass}_vs_sofue2020.png"

    echo
    echo "running a0=${a0}, M_HORB=${horb_mass}"

    cargo run -q -p curve_fitter -- total 3d_z2 "$a0" "$horb_mass" 6e10 3.0 1e10 0.7 > "$csv"

    ./scripts/plot_total_vs_mw.py \
      "$csv" \
      "$MW" \
      -o "$png" \
      --title "HORB 3d_z2 a0=${a0}, M=${horb_mass} + toy baryons vs Sofue 2020"

    ./scripts/score_total_vs_mw.py \
      "$csv" \
      "$MW" \
      --r-min 5 \
      --r-max 25 \
      | tail -n 1 >> "$REPORT"
  done
done

echo
echo "scores:"
column -s, -t "$REPORT" || cat "$REPORT"

echo
echo "plots:"
ls -lh plots/dz2_total_a*_vs_sofue2020.png
