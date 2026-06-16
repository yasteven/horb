#!/usr/bin/env python3

import argparse
from pathlib import Path
import re
import numpy as np


def extract_a0(path):
    m = re.search(r"_a([0-9.]+)_m", str(path))
    if not m:
        return np.nan
    return float(m.group(1))


def main():
    parser = argparse.ArgumentParser(description="Rank HORB x-z lobe extent against target Fermi-bubble height.")
    parser.add_argument("csv", type=Path)
    parser.add_argument("--target-total-height", type=float, default=18.0)
    parser.add_argument("--threshold", type=float, default=0.15)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True, dtype=None, encoding=None)

    rows = []
    for row in data:
        if abs(float(row["threshold"]) - args.threshold) > 1e-9:
            continue

        height = float(row["height_total_kpc"])
        if not np.isfinite(height):
            continue

        csv_path = str(row["csv"])
        a0 = extract_a0(csv_path)
        score = abs(height - args.target_total_height)

        rows.append((score, a0, height, float(row["z_top_kpc"]), float(row["z_bottom_kpc"]), csv_path))

    rows.sort(key=lambda x: x[0])

    print("rank,a0_kpc,threshold,total_height_kpc,z_top_kpc,z_bottom_kpc,abs_error_kpc,csv")
    for i, (score, a0, height, ztop, zbot, csv_path) in enumerate(rows, start=1):
        print(f"{i},{a0:.6f},{args.threshold:.6f},{height:.6f},{ztop:.6f},{zbot:.6f},{score:.6f},{csv_path}")


if __name__ == "__main__":
    main()
