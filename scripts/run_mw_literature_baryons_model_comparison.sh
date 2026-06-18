#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

./scripts/make_mw_sofue2020_curve.py
./scripts/convert_sofue2020_to_standard_rc.py

A0="${1:-1.0}"
HORB_MASS="${2:-2.5e11}"
BASELINE_MASS="${3:-2.5e11}"

RC="data/milky_way/sofue_2020_standard_rc.csv"
BARYONS="data/milky_way/sofue_literature_baryons.csv"

DM_COMPARE_CSV="curves/mw_lit_baryons_dm_compare_a${A0}_m${HORB_MASS}_base${BASELINE_MASS}.csv"
PLOT="plots/mw_literature_baryons_total_models_vs_sofue2020.png"

cargo run -q -p curve_fitter -- compare_dm_fixed_baselines \
  3d_z2 "$A0" "$HORB_MASS" "$BASELINE_MASS" \
  5.0 15.0 8.0 80.0 \
  > "$DM_COMPARE_CSV"

./scripts/plot_total_model_comparison_vs_standard_rc.py \
  --dm-compare-csv "$DM_COMPARE_CSV" \
  --rc-csv "$RC" \
  --baryons-csv "$BARYONS" \
  -o "$PLOT" \
  --title "Milky Way: HORB vs classical halos + literature baryons"

cargo run -q -p curve_fitter -- fit_standard_rc_multi_baryons \
  "$RC" "$BARYONS" \
  3d_z2 "$A0" "$HORB_MASS" \
  5 25 \
  > reports/mw_literature_baryons_horb_score_inner.csv

cargo run -q -p curve_fitter -- fit_standard_rc_multi_baryons \
  "$RC" "$BARYONS" \
  3d_z2 "$A0" "$HORB_MASS" \
  5 95.56 \
  > reports/mw_literature_baryons_horb_score_all.csv

echo
echo "wrote plot:"
ls -lh "$PLOT"

echo
echo "inner score:"
cat reports/mw_literature_baryons_horb_score_inner.csv

echo
echo "all score:"
cat reports/mw_literature_baryons_horb_score_all.csv
