#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Plot HORB vs classical DM halo comparison.")
    parser.add_argument("csv", type=Path)
    parser.add_argument("-o", "--out", type=Path, default=None)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True)

    r = data["r_kpc"]

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(r, data["v_horb_kms"], label="HORB 3d_z2")
    ax.plot(r, data["v_piso_kms"], label="pseudo-isothermal")
    ax.plot(r, data["v_nfw_kms"], label="NFW")
    ax.plot(r, data["v_burkert_kms"], label="Burkert")

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v_DM [km/s]")

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
