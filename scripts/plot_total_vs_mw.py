#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Plot HORB total curve against Milky Way rotation data.")
    parser.add_argument("model_csv", type=Path)
    parser.add_argument("mw_csv", type=Path)
    parser.add_argument("-o", "--out", type=Path, default=None)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    model = np.genfromtxt(args.model_csv, delimiter=",", names=True)
    mw = np.genfromtxt(args.mw_csv, delimiter=",", names=True, dtype=None, encoding=None)

    fig, ax = plt.subplots(figsize=(9, 6))

    r = model["r_kpc"]

    ax.plot(r, model["v_dm_kms"], label="HORB DM")
    ax.plot(r, model["v_disk_kms"], label="toy disk")
    ax.plot(r, model["v_bulge_kms"], label="toy bulge")
    ax.plot(r, model["v_total_kms"], label="HORB + toy baryons", linewidth=2.5)

    ax.errorbar(
        mw["R_kpc"],
        mw["v_kms"],
        yerr=mw["v_err_kms"],
        fmt="o",
        capsize=3,
        label="Milky Way target",
    )

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v [km/s]")

    if args.title is None:
        args.title = "HORB total curve vs Milky Way target"

    ax.set_title(args.title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if args.out is None:
        args.out = Path("plots") / f"{args.model_csv.stem}_vs_mw.png"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
