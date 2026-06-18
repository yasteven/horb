#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way reports plots curves

cd cuda_kernels/cuda
./build.sh
cd ../..

./scripts/make_mw_sofue2020_curve.py
./scripts/convert_sofue2020_to_standard_rc.py

RC="${1:-data/milky_way/sofue_2020_standard_rc.csv}"
BARYONS="${2:-data/milky_way/sofue_literature_baryons.csv}"
R_MIN="${3:-5}"
R_MAX="${4:-25}"
N_SIDE="${5:-64}"
EXTENT="${6:-80}"

OUT="reports/cuda_mw_wavefunction_scan_n${N_SIDE}_r${R_MIN}_${R_MAX}.csv"

HORB_WAVE_STATE_LIST="${HORB_WAVE_STATE_LIST:-3d_z2,3d_xy,3d_xz,3d_yz}" \
HORB_A0_LIST="${HORB_A0_LIST:-0.8,0.9,1.0,1.1,1.2,1.4}" \
HORB_DM_MASS_LIST="${HORB_DM_MASS_LIST:-1.5e11,2e11,2.5e11,3e11}" \
HORB_SOFTENING_LIST="${HORB_SOFTENING_LIST:-0.5}" \
HORB_COEFF_LEVELS="${HORB_COEFF_LEVELS:-2}" \
LD_LIBRARY_PATH="$PWD/cuda_kernels/cuda:${LD_LIBRARY_PATH:-}" \
cargo run -p cuda_kernels --bin scan_horb_cuda_mw_wavefunction -- \
  "$RC" \
  "$BARYONS" \
  "$OUT" \
  "$R_MIN" \
  "$R_MAX" \
  "$N_SIDE" \
  "$EXTENT"

echo
echo "best by chi2:"
(head -n 1 "$OUT" && tail -n +2 "$OUT" | sort -t, -k13,13n | head -30) | column -s, -t

echo
echo "best by RMS:"
(head -n 1 "$OUT" && tail -n +2 "$OUT" | sort -t, -k11,11n | head -30) | column -s, -t

echo
echo "wrote $OUT"
