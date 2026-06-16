#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Plot Milky Way rotation curve CSV.")
    parser.add_argument("csv", type=Path)
    parser.add_argument("-o", "--out", type=Path, default=None)
    parser.add_argument("--title", default="Milky Way rotation curve")
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True, dtype=None, encoding=None)

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.errorbar(
        data["R_kpc"],
        data["v_kms"],
        yerr=data["v_err_kms"],
        fmt="o-",
        capsize=3,
        linewidth=1.5,
        markersize=3,
        label="MW RC",
    )

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v [km/s]")
    ax.set_title(args.title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()

    if args.out is None:
        args.out = Path("plots") / f"{args.csv.stem}.png"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
