#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/milky_way reports plots images curves

cd cuda_kernels/cuda
./build.sh
cd ../..

./scripts/make_mw_sofue2020_curve.py
./scripts/convert_sofue2020_to_standard_rc.py

RC="data/milky_way/sofue_2020_standard_rc.csv"
BARYONS="data/milky_way/sofue_literature_baryons.csv"

echo
echo "=== 1. Build CUDA orbital basis library, n64 ==="

HORB_STATE_LIST="${HORB_STATE_LIST:-1s,2p_z,2p_x,2p_y,3d_z2,3d_x2_y2,3d_xy,3d_xz,3d_yz}" \
HORB_A0_LIST="${HORB_A0_LIST:-0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4}" \
LD_LIBRARY_PATH="$PWD/cuda_kernels/cuda:${LD_LIBRARY_PATH:-}" \
cargo run -p cuda_kernels --bin write_cuda_orbital_basis_library -- \
  reports/cuda_orbital_basis_library_n64.csv \
  64 \
  80 \
  0.5 \
  1e11

echo
echo "=== 2. Sparse positive density-basis fit ==="

./scripts/fit_cuda_orbital_basis_wavelets.py \
  --basis-csv reports/cuda_orbital_basis_library_n64.csv \
  --rc-csv "$RC" \
  --baryons-csv "$BARYONS" \
  --r-min 5 \
  --r-max 25 \
  --max-components 4 \
  --top-basis 24 \
  --metric chi2 \
  -o plots/cuda_orbital_basis_wavelet_fit_n64.png \
  --summary-out reports/cuda_orbital_basis_wavelet_fit_n64.csv

echo
echo "=== 3. True real wavefunction coefficient scan ==="

HORB_WAVE_STATE_LIST="${HORB_WAVE_STATE_LIST:-3d_z2,3d_xy,3d_xz,3d_yz}" \
HORB_A0_LIST="${HORB_WAVE_A0_LIST:-0.8,0.9,1.0,1.1,1.2,1.4}" \
HORB_DM_MASS_LIST="${HORB_WAVE_DM_MASS_LIST:-1.5e11,2e11,2.5e11,3e11}" \
HORB_SOFTENING_LIST="${HORB_WAVE_SOFTENING_LIST:-0.5}" \
HORB_COEFF_LEVELS="${HORB_COEFF_LEVELS:-2}" \
./scripts/scan_cuda_mw_wavefunction.sh \
  "$RC" \
  "$BARYONS" \
  5 25 64 80

echo
echo "=== 4. Wavefunction scan plots ==="

./scripts/plot_cuda_wavefunction_scan.py \
  --scan-csv reports/cuda_mw_wavefunction_scan_n64_r5_25.csv \
  --rc-csv "$RC" \
  --baryons-csv "$BARYONS" \
  --summary-out reports/cuda_mw_wavefunction_scan_n64_summary.csv \
  --bar-out plots/cuda_mw_wavefunction_scan_n64_top_chi2.png \
  --curve-out plots/cuda_mw_wavefunction_scan_n64_best_info.png \
  --metric chi2_per_point \
  --top 20

echo
echo "=== 5. Copy images ==="

cp plots/cuda_orbital_basis_wavelet_fit_n64.png \
   images/cuda_orbital_basis_wavelet_fit_n64.png

cp plots/cuda_mw_wavefunction_scan_n64_top_chi2.png \
   images/cuda_mw_wavefunction_scan_n64_top_chi2.png

cp plots/cuda_mw_wavefunction_scan_n64_best_info.png \
   images/cuda_mw_wavefunction_scan_n64_best_info.png

echo
echo "=== best sparse density-basis rows ==="
column -s, -t reports/cuda_orbital_basis_wavelet_fit_n64.csv | head -20 || true

echo
echo "=== best real-wavefunction rows ==="
column -s, -t reports/cuda_mw_wavefunction_scan_n64_summary.csv | head -20 || true

echo
echo "wrote:"
ls -lh \
  plots/cuda_orbital_basis_wavelet_fit_n64.png \
  plots/cuda_mw_wavefunction_scan_n64_top_chi2.png \
  plots/cuda_mw_wavefunction_scan_n64_best_info.png
