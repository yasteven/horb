#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Plot HORB total rotation curve CSV.")
    parser.add_argument("csv", type=Path)
    parser.add_argument("-o", "--out", type=Path, default=None)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True)

    r = data["r_kpc"]
    v_dm = data["v_dm_kms"]
    v_disk = data["v_disk_kms"]
    v_bulge = data["v_bulge_kms"]
    v_baryon = data["v_baryon_kms"]
    v_total = data["v_total_kms"]

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(r, v_dm, label="HORB DM")
    ax.plot(r, v_disk, label="disk")
    ax.plot(r, v_bulge, label="bulge")
    ax.plot(r, v_baryon, label="baryons")
    ax.plot(r, v_total, label="total", linewidth=2.5)

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v [km/s]")

    if args.title is None:
        args.title = args.csv.stem

    ax.set_title(args.title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if args.out is None:
        args.out = Path("plots") / f"{args.csv.stem}.png"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
