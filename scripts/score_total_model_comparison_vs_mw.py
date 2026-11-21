#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np


from horb_common import has_einasto, score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--total-csv", required=True, type=Path)
    parser.add_argument("--dm-compare-csv", required=True, type=Path)
    parser.add_argument("--mw-csv", required=True, type=Path)
    parser.add_argument("--r-min", type=float, default=5.0)
    parser.add_argument("--r-max", type=float, default=95.56)
    args = parser.parse_args()

    total = np.genfromtxt(args.total_csv, delimiter=",", names=True)
    dm = np.genfromtxt(args.dm_compare_csv, delimiter=",", names=True)
    mw = np.genfromtxt(args.mw_csv, delimiter=",", names=True, dtype=None, encoding=None)

    r = total["r_kpc"]
    v_baryon2 = total["v_disk_kms"]**2 + total["v_bulge_kms"]**2
    common = (r >= dm["r_kpc"][0]) & (r <= dm["r_kpc"][-1])
    r, v_baryon2 = r[common], v_baryon2[common]

    models = {
        "HORB+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_horb_kms"])**2 + v_baryon2),
        "NFW+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_nfw_kms"])**2 + v_baryon2),
        "pISO+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_piso_kms"])**2 + v_baryon2),
        "Burkert+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_burkert_kms"])**2 + v_baryon2),
        "baryons_only": np.sqrt(v_baryon2),
    }
    if has_einasto(dm):
        models["Einasto+baryons"] = np.sqrt(np.interp(r, dm["r_kpc"], dm["v_einasto_kms"])**2 + v_baryon2)

    print("model,r_min,r_max,n,rms_kms,chi2,chi2_per_point")
    for name, v in models.items():
        row = score(name, r, v, mw, args.r_min, args.r_max)
        print(f"{row[0]},{row[1]:.3f},{row[2]:.3f},{row[3]},{row[4]:.6f},{row[5]:.6f},{row[6]:.6f}")


if __name__ == "__main__":
    main()
