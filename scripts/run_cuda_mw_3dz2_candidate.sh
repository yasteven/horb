#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way curves plots reports

cd cuda_kernels/cuda
./build.sh
cd ../..

./scripts/make_mw_sofue2020_curve.py
./scripts/convert_sofue2020_to_standard_rc.py

A0="${1:-1.0}"
MASS="${2:-2.5e11}"
N_SIDE="${3:-96}"
EXTENT="${4:-80}"
SOFTENING="${5:-0.25}"

CUDA_CSV="curves/cuda_horb_diskplane_a$(printf "%.3f" "$A0")_m$(printf "%.3e" "$MASS")_n${N_SIDE}.csv"
SPH_CSV="curves/spherical_horb_compare_a${A0}_m${MASS}.csv"
PLOT="plots/mw_cuda_horb_diskplane_total_vs_sofue2020_a${A0}_m${MASS}_n${N_SIDE}.png"

LD_LIBRARY_PATH="$PWD/cuda_kernels/cuda:${LD_LIBRARY_PATH:-}" \
cargo run -p cuda_kernels --bin test_horb_disk_plane_curve -- \
  "$A0" "$MASS" "$N_SIDE" "$EXTENT" "$SOFTENING"

cargo run -q -p curve_fitter -- compare_dm_fixed_baselines \
  3d_z2 "$A0" "$MASS" "$MASS" \
  5.0 15.0 8.0 80.0 \
  > "$SPH_CSV"

./scripts/plot_cuda_total_vs_mw.py \
  --cuda-csv "$CUDA_CSV" \
  --dm-compare-csv "$SPH_CSV" \
  --rc-csv data/milky_way/sofue_2020_standard_rc.csv \
  --baryons-csv data/milky_way/sofue_literature_baryons.csv \
  -o "$PLOT" \
  --title "Milky Way CUDA HORB 3d_z2: a0=${A0}, M=${MASS}, n=${N_SIDE}"

echo
echo "wrote plot:"
ls -lh "$PLOT"
