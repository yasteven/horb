#!/usr/bin/env bash
set -euo pipefail

mkdir -p reports

RC="${1:-data/milky_way/sofue_2020_standard_rc.csv}"
BARYONS="${2:-data/milky_way/sofue_literature_baryons.csv}"
R_MIN="${3:-5}"
R_MAX="${4:-25}"

OUT="reports/mw_multi_baryons_horb_scan.csv"

echo "rc_csv,baryons_csv,state,a0_star_kpc,dm_mass_msun,r_min_kpc,r_max_kpc,n,rms_kms,chi2,chi2_per_point,has_rc_baryons,n_baryon_components" > "$OUT"

for a0 in 1.0 1.1 1.2 1.3 1.4 1.5 1.6; do
  for dm_mass in 2e11 2.5e11 3e11 3.5e11 4e11 4.5e11 5e11; do
    cargo run -q -p curve_fitter -- fit_standard_rc_multi_baryons \
      "$RC" \
      "$BARYONS" \
      3d_z2 "$a0" "$dm_mass" \
      "$R_MIN" "$R_MAX" \
      | tail -n 1 >> "$OUT"
  done
done

echo
echo "best by RMS:"
(head -n 1 "$OUT" && tail -n +2 "$OUT" | sort -t, -k9,9n | head -20) | column -s, -t

echo
echo
echo "best by chi2:"
(head -n 1 "$OUT" && tail -n +2 "$OUT" | sort -t, -k11,11n | head -20) | column -s, -t
