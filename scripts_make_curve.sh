#!/usr/bin/env bash
set -euo pipefail

cargo fmt
cargo test -p orbital_basis
cargo run -p curve_fitter > dz2_curve.csv

echo
echo "wrote dz2_curve.csv"

echo
echo "===== first rows ====="
head dz2_curve.csv

echo
echo "===== last rows ====="
tail dz2_curve.csv

echo
echo "===== max velocity row ====="
awk -F, '
NR == 1 { next }
$3 > vmax {
  vmax = $3
  r = $1
  m = $2
}
END {
  printf "r_kpc=%s  M_enc_Msun=%s  v_circ_kms=%.6f\n", r, m, vmax
}
' dz2_curve.csv
