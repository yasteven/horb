#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np


def summarize(name, r, v):
    i = np.nanargmax(v)
    return {
        "name": name,
        "vmax": v[i],
        "r_vmax": r[i],
        "v8": np.interp(8.0, r, v),
        "v20": np.interp(20.0, r, v),
        "v50": np.interp(50.0, r, v),
        "v80": np.interp(80.0, r, v),
    }


def main():
    parser = argparse.ArgumentParser(description="Summarize HORB vs classical halo comparison.")
    parser.add_argument("csv", type=Path)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True)
    r = data["r_kpc"]

    rows = [
        summarize("HORB", r, data["v_horb_kms"]),
        summarize("pseudo-isothermal", r, data["v_piso_kms"]),
        summarize("NFW", r, data["v_nfw_kms"]),
        summarize("Burkert", r, data["v_burkert_kms"]),
    ]

    print("model,vmax_kms,r_vmax_kpc,v8_kms,v20_kms,v50_kms,v80_kms,decline_20_to_80_kms")
    for row in rows:
        decline = row["v20"] - row["v80"]
        print(
            f"{row['name']},"
            f"{row['vmax']:.6f},"
            f"{row['r_vmax']:.6f},"
            f"{row['v8']:.6f},"
            f"{row['v20']:.6f},"
            f"{row['v50']:.6f},"
            f"{row['v80']:.6f},"
            f"{decline:.6f}"
        )


if __name__ == "__main__":
    main()
