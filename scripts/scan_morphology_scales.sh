#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves plots reports

MASS="${1:-1e11}"
REPORT="reports/morphology_scale_scan_m${MASS}.csv"

echo "a0_kpc,density_peak_kpc,curve_vmax_radius_kpc,curve_menc_at_vmax_msun,curve_vmax_kms,rho_max_x_kpc,rho_max_z_kpc,rho_max_r_kpc,rho_max_msun_per_kpc3" > "$REPORT"

for a0 in 1.0 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.0; do
  echo
  echo "===== a0=${a0} kpc, M=${MASS} M_sun ====="

  CURVE="curves/dz2_curve_a${a0}_m${MASS}.csv"
  XZ="curves/dz2_xz_a${a0}_m${MASS}.csv"

  cargo run -q -p curve_fitter -- curve 3d_z2 "$a0" "$MASS" > "$CURVE"
  cargo run -q -p curve_fitter -- xz 3d_z2 "$a0" "$MASS" > "$XZ"

  ./scripts/plot_xz_density.py "$XZ" \
    --mode log \
    --contours \
    -o "plots/dz2_xz_a${a0}_m${MASS}_log.png" \
    --title "3d_z2 log density, a0=${a0} kpc, M=${MASS} M_sun"

  ./scripts/plot_xz_density.py "$XZ" \
    --mode linear \
    --contours \
    -o "plots/dz2_xz_a${a0}_m${MASS}_linear.png" \
    --title "3d_z2 normalized density, a0=${a0} kpc, M=${MASS} M_sun"

  ./scripts/plot_xz_density.py "$XZ" \
    --mode threshold \
    --threshold 0.15 \
    --contours \
    -o "plots/dz2_xz_a${a0}_m${MASS}_threshold015.png" \
    --title "3d_z2 threshold density > 0.15 max, a0=${a0} kpc"

  CURVE_SUMMARY="$(awk -F, '
    NR == 1 { next }
    $3 > vmax {
      vmax = $3
      r_vmax = $1
      m_vmax = $2
    }
    END {
      printf "%s,%s,%.6f", r_vmax, m_vmax, vmax
    }
  ' "$CURVE")"

  XZ_SUMMARY="$(python3 - <<PY
import numpy as np
data = np.genfromtxt("$XZ", delimiter=",", names=True)
rho = data["rho_Msun_per_kpc3"]
i = np.nanargmax(rho)
print(f'{data["x_kpc"][i]:.6f},{data["z_kpc"][i]:.6f},{data["r_kpc"][i]:.6f},{rho[i]:.6e}')
PY
)"

  DENSITY_PEAK="$(python3 - <<PY
a0 = float("$a0")
print(f"{6.0 * a0:.6f}")
PY
)"

  echo "${a0},${DENSITY_PEAK},${CURVE_SUMMARY},${XZ_SUMMARY}" >> "$REPORT"

  echo "density peak target: ${DENSITY_PEAK} kpc"
  echo "curve summary: ${CURVE_SUMMARY}"
  echo "xz max summary: ${XZ_SUMMARY}"
done

echo
echo "wrote report:"
echo "$REPORT"
column -s, -t "$REPORT" || cat "$REPORT"
