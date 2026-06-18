#!/usr/bin/env bash
set -euo pipefail

mkdir -p reports

# Override these after first scan if needed.
export HORB_A0_LIST="${HORB_A0_LIST:-0.9,1.0,1.1}"
export HORB_DM_MASS_LIST="${HORB_DM_MASS_LIST:-2e11,2.5e11,3e11}"
export HORB_SOFTENING_LIST="${HORB_SOFTENING_LIST:-0.1,0.25,0.5}"

for n in 64 96 128; do
  echo
  echo "===== CUDA MW 3d_z2 convergence n_side=$n ====="
  ./scripts/scan_cuda_mw_3dz2.sh \
    data/milky_way/sofue_2020_standard_rc.csv \
    data/milky_way/sofue_literature_baryons.csv \
    5 \
    25 \
    "$n" \
    80
done
