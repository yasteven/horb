#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Measure HORB 3d_z2 x-z lobe extent from density CSV.")
    parser.add_argument("csv", type=Path)
    parser.add_argument("--thresholds", nargs="+", type=float, default=[0.05, 0.10, 0.15, 0.20, 0.25, 0.50])
    parser.add_argument("--min-polar-z", type=float, default=0.0)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True)

    x = data["x_kpc"]
    z = data["z_kpc"]
    r = data["r_kpc"]
    rho = data["rho_Msun_per_kpc3"]

    rho_max = np.nanmax(rho)
    norm = rho / rho_max

    print("csv,threshold,z_top_kpc,z_bottom_kpc,height_total_kpc,x_at_top_kpc,x_at_bottom_kpc,r_at_top_kpc,r_at_bottom_kpc")

    for th in args.thresholds:
        mask_top = (norm >= th) & (z >= args.min_polar_z)
        mask_bottom = (norm >= th) & (z <= -args.min_polar_z)

        if not np.any(mask_top) or not np.any(mask_bottom):
            print(f"{args.csv},{th:.6f},nan,nan,nan,nan,nan,nan,nan")
            continue

        i_top_candidates = np.where(mask_top)[0]
        i_bottom_candidates = np.where(mask_bottom)[0]

        i_top = i_top_candidates[np.argmax(z[i_top_candidates])]
        i_bottom = i_bottom_candidates[np.argmin(z[i_bottom_candidates])]

        z_top = z[i_top]
        z_bottom = z[i_bottom]
        height = z_top - z_bottom

        print(
            f"{args.csv},"
            f"{th:.6f},"
            f"{z_top:.6f},"
            f"{z_bottom:.6f},"
            f"{height:.6f},"
            f"{x[i_top]:.6f},"
            f"{x[i_bottom]:.6f},"
            f"{r[i_top]:.6f},"
            f"{r[i_bottom]:.6f}"
        )


if __name__ == "__main__":
    main()
