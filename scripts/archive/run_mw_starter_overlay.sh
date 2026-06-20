#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

./scripts/make_mw_starter_curve.py

echo "model_csv,r_min,r_max,n,rms_kms,chi2,chi2_per_point" > reports/score_horb_total_vs_mw_starter.csv

for a0 in 1.4 1.5 1.6 1.7; do
  csv="curves/dz2_total_a${a0}_m5e10_disk6e10_bulge1e10.csv"
  png="plots/dz2_total_a${a0}_m5e10_vs_mw_starter.png"

  cargo run -q -p curve_fitter -- total 3d_z2 "$a0" 5e10 6e10 3.0 1e10 0.7 > "$csv"

  ./scripts/plot_total_vs_mw.py \
    "$csv" \
    data/milky_way/mw_starter_total_curve.csv \
    -o "$png" \
    --title "HORB 3d_z2 a0=${a0}, M=5e10 + toy baryons vs MW starter"

  ./scripts/score_total_vs_mw.py \
    "$csv" \
    data/milky_way/mw_starter_total_curve.csv \
    --r-min 5 \
    --r-max 25 \
    | tail -n 1 >> reports/score_horb_total_vs_mw_starter.csv
done

echo
echo "scores:"
column -s, -t reports/score_horb_total_vs_mw_starter.csv || cat reports/score_horb_total_vs_mw_starter.csv

echo
echo "plots:"
ls -lh plots/dz2_total_a*_m5e10_vs_mw_starter.png
