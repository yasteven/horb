#!/usr/bin/env bash
set -euo pipefail

for a0 in 1.2 1.5 1.7 2.0; do
  ./scripts/plot_orbital_views.sh "$a0" 1e11
done

echo
echo "generated plots:"
ls -lh plots/dz2_xz_a*_m1e11_*.png
