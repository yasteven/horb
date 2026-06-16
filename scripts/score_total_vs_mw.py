#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Score total model curve against MW rotation target.")
    parser.add_argument("model_csv", type=Path)
    parser.add_argument("mw_csv", type=Path)
    parser.add_argument("--r-min", type=float, default=5.0)
    parser.add_argument("--r-max", type=float, default=25.0)
    args = parser.parse_args()

    model = np.genfromtxt(args.model_csv, delimiter=",", names=True)
    mw = np.genfromtxt(args.mw_csv, delimiter=",", names=True, dtype=None, encoding=None)

    mask = (mw["R_kpc"] >= args.r_min) & (mw["R_kpc"] <= args.r_max)

    r_obs = mw["R_kpc"][mask]
    v_obs = mw["v_kms"][mask]
    err = mw["v_err_kms"][mask]

    v_model = np.interp(r_obs, model["r_kpc"], model["v_total_kms"])

    resid = v_model - v_obs
    rms = np.sqrt(np.mean(resid**2))
    chi2 = np.sum((resid / err) ** 2)
    n = len(r_obs)

    print("model_csv,r_min,r_max,n,rms_kms,chi2,chi2_per_point")
    print(f"{args.model_csv},{args.r_min:.3f},{args.r_max:.3f},{n},{rms:.6f},{chi2:.6f},{chi2/n:.6f}")


if __name__ == "__main__":
    main()
