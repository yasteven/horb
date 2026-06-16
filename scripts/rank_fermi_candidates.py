#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Rank HORB morphology scan by target lobe height.")
    parser.add_argument("report", type=Path)
    parser.add_argument("--target-min", type=float, default=8.0)
    parser.add_argument("--target-max", type=float, default=10.5)
    args = parser.parse_args()

    data = np.genfromtxt(args.report, delimiter=",", names=True)

    target_mid = 0.5 * (args.target_min + args.target_max)
    peak = data["density_peak_kpc"]

    score = np.abs(peak - target_mid)
    order = np.argsort(score)

    print(f"target density lobe height: {args.target_min:.2f} to {args.target_max:.2f} kpc")
    print()
    print("rank,a0_kpc,density_peak_kpc,vmax_radius_kpc,vmax_kms,inside_target")

    for rank, i in enumerate(order, start=1):
        inside = args.target_min <= peak[i] <= args.target_max
        print(
            f"{rank},"
            f"{data['a0_kpc'][i]:.3f},"
            f"{data['density_peak_kpc'][i]:.3f},"
            f"{data['curve_vmax_radius_kpc'][i]:.3f},"
            f"{data['curve_vmax_kms'][i]:.3f},"
            f"{inside}"
        )


if __name__ == "__main__":
    main()
