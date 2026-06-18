#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np


def score(name, r_model, v_model, mw, r_min, r_max):
    mask = (mw["R_kpc"] >= r_min) & (mw["R_kpc"] <= r_max)
    r_obs = mw["R_kpc"][mask]

    # Accept either old Sofue helper schema or standard RC schema.
    if "Vobs_kms" in mw.dtype.names:
        v_obs = mw["Vobs_kms"][mask]
        err = mw["eV_kms"][mask]
    else:
        v_obs = mw["v_kms"][mask]
        err = mw["v_err_kms"][mask]

    v_interp = np.interp(r_obs, r_model, v_model)
    resid = v_interp - v_obs
    rms = np.sqrt(np.mean(resid**2))
    chi2 = np.sum((resid / err) ** 2)
    n = len(r_obs)

    return name, r_min, r_max, n, rms, chi2, chi2 / n


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

    models = {
        "HORB+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_horb_kms"])**2 + v_baryon2),
        "NFW+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_nfw_kms"])**2 + v_baryon2),
        "pISO+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_piso_kms"])**2 + v_baryon2),
        "Burkert+baryons": np.sqrt(np.interp(r, dm["r_kpc"], dm["v_burkert_kms"])**2 + v_baryon2),
        "baryons_only": np.sqrt(v_baryon2),
    }

    print("model,r_min,r_max,n,rms_kms,chi2,chi2_per_point")
    for name, v in models.items():
        row = score(name, r, v, mw, args.r_min, args.r_max)
        print(f"{row[0]},{row[1]:.3f},{row[2]:.3f},{row[3]},{row[4]:.6f},{row[5]:.6f},{row[6]:.6f}")


if __name__ == "__main__":
    main()
