#!/usr/bin/env bash
set -euo pipefail

mkdir -p curves

for a0 in 0.8 1.0 1.2 1.5 1.7; do
  out="curves/dz2_density_axes_a${a0}_m1e11.csv"

  echo "running density axes a0_star=${a0} kpc -> ${out}"
  cargo run -q -p curve_fitter -- density 3d_z2 "$a0" 1e11 > "$out"

  awk -F, -v a0="$a0" '
  NR == 1 { next }
  $2 > max_pole {
    max_pole = $2
    r_pole = $1
  }
  $3 > max_eq {
    max_eq = $3
    r_eq = $1
  }
  END {
    printf "a0=%s kpc  expected_radial_peak≈%.3f kpc  pole_peak_r=%s kpc  equator_peak_r=%s kpc  pole/equator=%.3f\n",
      a0, 9*a0, r_pole, r_eq, max_pole/max_eq
  }
  ' "$out"
done
