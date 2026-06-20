#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

A0="1.5"
HORB_MASS="2.5e11"
DISK_MASS="1e11"
DISK_SCALE="3.0"
BULGE_MASS="1e10"
BULGE_SCALE="0.7"

./scripts/make_mw_sofue2020_curve.py

MW="data/milky_way/sofue_2020_unified_rc.csv"

TOTAL_CSV="curves/best_candidate_total.csv"
TOTAL_PNG="plots/best_candidate_total_vs_sofue2020.png"

XZ_CSV="curves/best_candidate_xz_ext40.csv"

echo "Generating best candidate total curve..."
cargo run -q -p curve_fitter -- total 3d_z2 "$A0" "$HORB_MASS" "$DISK_MASS" "$DISK_SCALE" "$BULGE_MASS" "$BULGE_SCALE" > "$TOTAL_CSV"

echo "Plotting best candidate vs Sofue 2020..."
./scripts/plot_total_vs_mw.py \
  "$TOTAL_CSV" \
  "$MW" \
  -o "$TOTAL_PNG" \
  --title "Best candidate: HORB 3d_z2 a0=${A0}, M=${HORB_MASS}, disk=${DISK_MASS}"

echo "Scoring best candidate..."
./scripts/score_total_vs_mw.py "$TOTAL_CSV" "$MW" --r-min 5 --r-max 25 \
  > reports/best_candidate_score_inner.csv

./scripts/score_total_vs_mw.py "$TOTAL_CSV" "$MW" --r-min 5 --r-max 95.56 \
  > reports/best_candidate_score_all.csv

echo "Generating morphology slice..."
HORB_XZ_EXTENT_KPC=40 HORB_XZ_N=321 \
  cargo run -q -p curve_fitter -- xz 3d_z2 "$A0" "$HORB_MASS" > "$XZ_CSV"

for TH in 0.15 0.20 0.25; do
  TAG="$(printf "%.2f" "$TH" | tr -d '.')"
  ./scripts/plot_xz_density.py "$XZ_CSV" \
    --mode threshold \
    --threshold "$TH" \
    --contours \
    -o "plots/best_candidate_xz_threshold${TAG}.png" \
    --title "Best candidate x-z morphology, threshold ${TH}"
done

./scripts/measure_xz_lobe_extent.py "$XZ_CSV" \
  --thresholds 0.10 0.15 0.20 0.25 0.30 0.40 0.50 \
  > reports/best_candidate_lobe_extent.csv

echo
echo "Best candidate generated:"
ls -lh "$TOTAL_PNG" plots/best_candidate_xz_threshold*.png

echo
echo "Inner score:"
cat reports/best_candidate_score_inner.csv

echo
echo "All-range score:"
cat reports/best_candidate_score_all.csv

echo
echo "Lobe extent:"
cat reports/best_candidate_lobe_extent.csv
