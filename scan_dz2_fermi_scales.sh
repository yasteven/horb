#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves

for a0 in 0.8 1.0 1.2 1.5 1.7 2.5 8.5; do
  out="curves/dz2_a${a0}_m1e11.csv"

  echo "running a0_star=${a0} kpc -> ${out}"
  cargo run -q -p curve_fitter -- 3d_z2 "$a0" 1e11 > "$out"

  awk -F, -v a0="$a0" '
  NR == 1 { next }
  $3 > vmax {
    vmax = $3
    r_vmax = $1
    m_vmax = $2
  }
  END {
    printf "a0=%s kpc  radial_density_peak≈%.3f kpc  vmax_radius=%s kpc  M_enc=%s  vmax=%.3f km/s\n",
      a0, 9*a0, r_vmax, m_vmax, vmax
  }
  ' "$out"
done
