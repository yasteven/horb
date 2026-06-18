#!/usr/bin/env bash
set -euo pipefail

RC="${1:-data/milky_way/sofue_2020_standard_rc.csv}"
R_MIN="${2:-5}"
R_MAX="${3:-25}"

mkdir -p reports

OUT="reports/standard_rc_scan_horb.csv"
echo "csv,state,a0_star_kpc,dm_mass_msun,disk_mass_msun,disk_scale_kpc,bulge_mass_msun,bulge_scale_kpc,r_min_kpc,r_max_kpc,n,rms_kms,chi2,chi2_per_point,has_baryons" > "$OUT"

for a0 in 1.3 1.4 1.5 1.6 1.7; do
  for dm_mass in 2e11 2.5e11 3e11; do
    for disk_mass in 8e10 1e11 1.2e11; do
      cargo run -q -p curve_fitter -- fit_standard_rc \
        "$RC" \
        3d_z2 "$a0" "$dm_mass" \
        "$disk_mass" 3.0 \
        1e10 0.7 \
        "$R_MIN" "$R_MAX" \
        | tail -n 1 >> "$OUT"
    done
  done
done

echo
echo "best fits:"
(head -n 1 "$OUT" && tail -n +2 "$OUT" | sort -t, -k12,12n | head -20) | column -s, -t
