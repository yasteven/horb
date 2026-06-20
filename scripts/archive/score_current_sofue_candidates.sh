#!/usr/bin/env bash
set -euo pipefail

MW="data/milky_way/sofue_2020_unified_rc.csv"
OUT="reports/score_current_sofue_candidates_ranges.csv"

echo "model_csv,range_name,r_min,r_max,n,rms_kms,chi2,chi2_per_point" > "$OUT"

for csv in curves/dz2_total_a1.5_m*_disk*_bulge1e10.csv; do
  [ -f "$csv" ] || continue

  for spec in "inner,5,25" "mid,10,40" "outer,25,95.56" "all,5,95.56"; do
    IFS=, read -r name rmin rmax <<< "$spec"

    line="$(
      ./scripts/score_total_vs_mw.py "$csv" "$MW" --r-min "$rmin" --r-max "$rmax" \
        | tail -n 1
    )"

    echo "$line" | awk -F, -v name="$name" '
      {
        printf "%s,%s,%s,%s,%s,%s,%s,%s\n", $1, name, $2, $3, $4, $5, $6, $7
      }
    ' >> "$OUT"
  done
done

echo
echo "best inner:"
awk -F, 'NR==1 || $2=="inner"' "$OUT" | sort -t, -k7,7n | head -10

echo
echo "best all:"
awk -F, 'NR==1 || $2=="all"' "$OUT" | sort -t, -k7,7n | head -10
