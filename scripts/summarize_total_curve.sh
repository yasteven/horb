#!/usr/bin/env bash
set -euo pipefail

CSV="$1"

echo "file: $CSV"

awk -F, '
NR == 1 { next }

$3 > vdm_max {
  vdm_max = $3
  r_vdm = $1
}

$4 > vdisk_max {
  vdisk_max = $4
  r_vdisk = $1
}

$5 > vbulge_max {
  vbulge_max = $5
  r_vbulge = $1
}

$6 > vbaryon_max {
  vbaryon_max = $6
  r_vbaryon = $1
}

$7 > vtotal_max {
  vtotal_max = $7
  r_vtotal = $1
}

END {
  printf "dm max:      r=%s kpc  v=%.3f km/s\n", r_vdm, vdm_max
  printf "disk max:    r=%s kpc  v=%.3f km/s\n", r_vdisk, vdisk_max
  printf "bulge max:   r=%s kpc  v=%.3f km/s\n", r_vbulge, vbulge_max
  printf "baryon max:  r=%s kpc  v=%.3f km/s\n", r_vbaryon, vbaryon_max
  printf "total max:   r=%s kpc  v=%.3f km/s\n", r_vtotal, vtotal_max
}
' "$CSV"
