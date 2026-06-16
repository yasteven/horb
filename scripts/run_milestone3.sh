#!/usr/bin/env bash
set -euo pipefail

MASS="${1:-1e11}"

cargo fmt
cargo test -p orbital_basis

./scripts/scan_morphology_scales.sh "$MASS"

./scripts/make_morphology_contact_sheet.py --mode threshold015 --mass "$MASS"
./scripts/make_morphology_contact_sheet.py --mode linear --mass "$MASS"
./scripts/make_morphology_contact_sheet.py --mode log --mass "$MASS"

./scripts/rank_fermi_candidates.py "reports/morphology_scale_scan_m${MASS}.csv" \
  --target-min 8.0 \
  --target-max 10.5

echo
echo "milestone 3 complete"
echo
echo "main report:"
echo "reports/morphology_scale_scan_m${MASS}.csv"
echo
echo "contact sheets:"
ls -lh plots/morphology_contact_sheet_*_m${MASS}.png
