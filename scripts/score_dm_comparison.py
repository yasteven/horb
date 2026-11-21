#!/usr/bin/env python3
"""Score every fixed halo with the same literature baryons and observations."""
import argparse
import csv
import sys
from pathlib import Path
import numpy as np
from horb_common import has_einasto, load_baryon_velocity_squared, score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dm-compare-csv", required=True, type=Path)
    parser.add_argument("--rc-csv", required=True, type=Path)
    parser.add_argument("--baryons-csv", required=True, type=Path)
    parser.add_argument("--cuda-csv", type=Path, help="also score the CUDA HORB curve on the common radius interval")
    parser.add_argument("--r-min", type=float, default=5.0)
    parser.add_argument("--r-max", type=float, default=80.0)
    args = parser.parse_args()
    dm = np.atleast_1d(np.genfromtxt(args.dm_compare_csv, delimiter=",", names=True))
    rc = np.atleast_1d(np.genfromtxt(args.rc_csv, delimiter=",", names=True, dtype=None, encoding=None))
    r = dm["r_kpc"]
    cuda = None
    if args.cuda_csv:
        cuda = np.atleast_1d(np.genfromtxt(args.cuda_csv, delimiter=",", names=True))
        lower, upper = max(r[0], cuda["r_kpc"][0]), min(r[-1], cuda["r_kpc"][-1])
        r = np.unique(np.concatenate((r, cuda["r_kpc"])))
        r = r[(r >= lower) & (r <= upper)]
    vb2 = load_baryon_velocity_squared(args.baryons_csv, r)
    models = [("HORB_spherical+baryons", "horb"), ("NFW+baryons", "nfw"),
              ("pISO+baryons", "piso"), ("Burkert+baryons", "burkert")]
    if has_einasto(dm):
        models.append(("Einasto+baryons", "einasto"))
    rows = [score(name, r, np.sqrt(np.interp(r, dm["r_kpc"], dm[f"v_{key}_kms"])**2 + vb2), rc, args.r_min, args.r_max)
            for name, key in models]
    if cuda is not None:
        rows.insert(0, score("HORB_CUDA+baryons", r,
                    np.sqrt(np.interp(r, cuda["r_kpc"], cuda["v_horb_cuda_kms"])**2 + vb2),
                    rc, args.r_min, args.r_max))
    rows.append(score("baryons_only", r, np.sqrt(vb2), rc, args.r_min, args.r_max))
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["model", "r_min", "r_max", "n", "rms_kms", "chi2", "chi2_per_point"])
    for name, lo, hi, n, rms, chi2, per_point in rows:
        writer.writerow([name, f"{lo:.3f}", f"{hi:.3f}", n, f"{rms:.6f}", f"{chi2:.6f}", f"{per_point:.6f}"])


if __name__ == "__main__":
    main()
