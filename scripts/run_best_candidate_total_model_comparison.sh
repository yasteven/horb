#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves plots reports data/milky_way

./scripts/make_mw_sofue2020_curve.py

A0="1.5"
HORB_MASS="2.5e11"
BASELINE_MASS="2.5e11"
DISK_MASS="1e11"
DISK_SCALE="3.0"
BULGE_MASS="1e10"
BULGE_SCALE="0.7"

MW="data/milky_way/sofue_2020_unified_rc.csv"

TOTAL_CSV="curves/best_candidate_total_a${A0}_m${HORB_MASS}_disk${DISK_MASS}.csv"
DM_COMPARE_CSV="curves/best_candidate_dm_compare_baselines.csv"
PLOT="plots/best_candidate_total_models_vs_sofue2020.png"

cargo run -q -p curve_fitter -- total 3d_z2 "$A0" "$HORB_MASS" "$DISK_MASS" "$DISK_SCALE" "$BULGE_MASS" "$BULGE_SCALE" \
  > "$TOTAL_CSV"

cargo run -q -p curve_fitter -- compare_dm_fixed_baselines 3d_z2 "$A0" "$HORB_MASS" "$BASELINE_MASS" 5.0 15.0 8.0 80.0 \
  > "$DM_COMPARE_CSV"

./scripts/plot_total_model_comparison_vs_mw.py \
  --total-csv "$TOTAL_CSV" \
  --dm-compare-csv "$DM_COMPARE_CSV" \
  --mw-csv "$MW" \
  -o "$PLOT" \
  --title "HORB vs classical halos + same baryons vs Sofue 2020"

./scripts/score_total_model_comparison_vs_mw.py \
  --total-csv "$TOTAL_CSV" \
  --dm-compare-csv "$DM_COMPARE_CSV" \
  --mw-csv "$MW" \
  --r-min 5 \
  --r-max 25 \
  > reports/best_candidate_total_models_score_inner.csv

./scripts/score_total_model_comparison_vs_mw.py \
  --total-csv "$TOTAL_CSV" \
  --dm-compare-csv "$DM_COMPARE_CSV" \
  --mw-csv "$MW" \
  --r-min 5 \
  --r-max 95.56 \
  > reports/best_candidate_total_models_score_all.csv

echo
echo "wrote plot:"
ls -lh "$PLOT"

echo
echo "inner score:"
cat reports/best_candidate_total_models_score_inner.csv

echo
echo "all score:"
cat reports/best_candidate_total_models_score_all.csv
