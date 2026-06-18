#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

./scripts/make_mw_sofue2020_curve.py
./scripts/convert_sofue2020_to_standard_rc.py

RC="data/milky_way/sofue_2020_standard_rc.csv"

# Sofue grand-RC bulge+disk decomposition values.
# M_b = 1.652e10 Msun, a_b = 0.522 kpc
# M_d = 3.41e10 Msun, a_d = 3.19 kpc
DISK_MASS="3.41e10"
DISK_SCALE="3.19"
BULGE_MASS="1.652e10"
BULGE_SCALE="0.522"

OUT="reports/mw_literature_baryons_horb_scan.csv"

echo "csv,state,a0_star_kpc,dm_mass_msun,disk_mass_msun,disk_scale_kpc,bulge_mass_msun,bulge_scale_kpc,r_min_kpc,r_max_kpc,n,rms_kms,chi2,chi2_per_point,has_baryons" > "$OUT"

for a0 in 1.2 1.3 1.4 1.5 1.6 1.7 1.8; do
  for dm_mass in 2e11 2.5e11 3e11 3.5e11 4e11 5e11 6e11 8e11; do
    cargo run -q -p curve_fitter -- fit_standard_rc \
      "$RC" \
      3d_z2 "$a0" "$dm_mass" \
      "$DISK_MASS" "$DISK_SCALE" \
      "$BULGE_MASS" "$BULGE_SCALE" \
      5 25 \
      | tail -n 1 >> "$OUT"
  done
done

echo
echo "best inner fits, Sofue literature baryons:"
(head -n 1 "$OUT" && tail -n +2 "$OUT" | sort -t, -k12,12n | head -20) | column -s, -t
