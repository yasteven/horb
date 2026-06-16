#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

echo
echo "============================================================"
echo "1. Generate Sofue 2020 Milky Way data"
echo "============================================================"
./scripts/make_mw_sofue2020_curve.py

echo
echo "============================================================"
echo "2. Plot Sofue 2020 Milky Way rotation curve"
echo "============================================================"
./scripts/plot_mw_curve.py \
  data/milky_way/sofue_2020_unified_rc.csv \
  -o plots/sofue_2020_unified_rc.png \
  --title "Sofue 2020 unified Milky Way rotation curve"

echo
echo "============================================================"
echo "3. Generate current best candidate"
echo "============================================================"

A0="1.5"
HORB_MASS="2.5e11"
DISK_MASS="1e11"
DISK_SCALE="3.0"
BULGE_MASS="1e10"
BULGE_SCALE="0.7"
MW="data/milky_way/sofue_2020_unified_rc.csv"

TOTAL_CSV="curves/best_candidate_total_a${A0}_m${HORB_MASS}_disk${DISK_MASS}.csv"
TOTAL_PNG="plots/best_candidate_total_vs_sofue2020.png"
XZ_CSV="curves/best_candidate_xz_a${A0}_m${HORB_MASS}_ext40.csv"

cargo run -q -p curve_fitter -- total 3d_z2 "$A0" "$HORB_MASS" "$DISK_MASS" "$DISK_SCALE" "$BULGE_MASS" "$BULGE_SCALE" \
  > "$TOTAL_CSV"

./scripts/plot_total_vs_mw.py \
  "$TOTAL_CSV" \
  "$MW" \
  -o "$TOTAL_PNG" \
  --title "Best candidate: HORB 3d_z2 a0=${A0}, M=${HORB_MASS}, disk=${DISK_MASS}"

echo
echo "============================================================"
echo "4. Score current best candidate"
echo "============================================================"

./scripts/score_total_vs_mw.py "$TOTAL_CSV" "$MW" --r-min 5 --r-max 25 \
  > reports/best_candidate_score_inner_5_25.csv

./scripts/score_total_vs_mw.py "$TOTAL_CSV" "$MW" --r-min 5 --r-max 95.56 \
  > reports/best_candidate_score_all_5_95.csv

echo "Inner score:"
cat reports/best_candidate_score_inner_5_25.csv

echo
echo "All score:"
cat reports/best_candidate_score_all_5_95.csv

echo
echo "============================================================"
echo "5. Generate best candidate Fermi morphology x-z slice"
echo "============================================================"

HORB_XZ_EXTENT_KPC=40 HORB_XZ_N=321 \
  cargo run -q -p curve_fitter -- xz 3d_z2 "$A0" "$HORB_MASS" \
  > "$XZ_CSV"

for TH in 0.15 0.20 0.25 0.50; do
  TAG="$(printf "%.2f" "$TH" | tr -d '.')"

  ./scripts/plot_xz_density.py "$XZ_CSV" \
    --mode threshold \
    --threshold "$TH" \
    --contours \
    -o "plots/best_candidate_xz_threshold${TAG}.png" \
    --title "Best candidate x-z morphology, rho/rho_max >= ${TH}"
done

./scripts/plot_xz_density.py "$XZ_CSV" \
  --mode log \
  --contours \
  -o "plots/best_candidate_xz_log.png" \
  --title "Best candidate x-z log density"

./scripts/plot_xz_density.py "$XZ_CSV" \
  --mode linear \
  --contours \
  -o "plots/best_candidate_xz_linear.png" \
  --title "Best candidate x-z normalized density"

./scripts/measure_xz_lobe_extent.py "$XZ_CSV" \
  --thresholds 0.10 0.15 0.20 0.25 0.30 0.40 0.50 \
  > reports/best_candidate_lobe_extent.csv

echo
echo "Best candidate lobe extent:"
cat reports/best_candidate_lobe_extent.csv

echo
echo "============================================================"
echo "6. Regenerate ext40 Fermi morphology scan/contact sheets"
echo "============================================================"

./scripts/scan_fermi_lobe_extent_ext40.sh "$HORB_MASS" 40 321

./scripts/make_fermi_ext40_contact_sheet.py --mass "$HORB_MASS" --extent 40 --threshold 015
./scripts/make_fermi_ext40_contact_sheet.py --mass "$HORB_MASS" --extent 40 --threshold 020
./scripts/make_fermi_ext40_contact_sheet.py --mass "$HORB_MASS" --extent 40 --threshold 025
./scripts/make_fermi_ext40_contact_sheet.py --mass "$HORB_MASS" --extent 40 --threshold 050

echo
echo "============================================================"
echo "7. Regenerate HORB/disk score table"
echo "============================================================"

./scripts/scan_horb_disk_vs_sofue2020.sh "$A0"
./scripts/score_current_sofue_candidates.sh

echo
echo "============================================================"
echo "8. Final output inventory"
echo "============================================================"

echo
echo "Core plots:"
ls -lh \
  plots/sofue_2020_unified_rc.png \
  plots/best_candidate_total_vs_sofue2020.png \
  plots/best_candidate_xz_threshold*.png \
  plots/best_candidate_xz_log.png \
  plots/best_candidate_xz_linear.png \
  plots/fermi_ext40_threshold*_m${HORB_MASS}_contact.png

echo
echo "Core reports:"
ls -lh \
  reports/best_candidate_score_inner_5_25.csv \
  reports/best_candidate_score_all_5_95.csv \
  reports/best_candidate_lobe_extent.csv \
  reports/fermi_lobe_extent_scan_m${HORB_MASS}_ext40.csv \
  reports/score_horb_disk_scan_a${A0}_vs_sofue2020.csv \
  reports/score_current_sofue_candidates_ranges.csv

echo
echo "Done."
