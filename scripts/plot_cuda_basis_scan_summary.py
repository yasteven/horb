#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--metric", choices=["rms_kms", "chi2_per_point"], default="chi2_per_point")
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True, dtype=None, encoding=None)

    states = sorted(set(str(s) for s in data["state"]))

    best = []
    for state in states:
        rows = data[data["state"] == state]
        idx = np.argmin(rows[args.metric])
        row = rows[idx]
        best.append((
            state,
            float(row[args.metric]),
            float(row["a0_star_kpc"]),
            float(row["dm_mass_msun"]),
            float(row["softening_kpc"]),
        ))

    best.sort(key=lambda x: x[1])

    labels = [x[0] for x in best]
    values = [x[1] for x in best]

    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ax.bar(labels, values)

    ax.set_xlabel("Hydrogenic orbital basis state")
    ax.set_ylabel(args.metric)
    ax.set_title(f"Best Milky Way CUDA disk-plane fit by orbital basis state ({args.metric})")
    ax.grid(axis="y", alpha=0.3)

    for i, (state, value, a0, mass, soft) in enumerate(best):
        ax.text(
            i,
            value,
            f"a0={a0:.2f}\nM={mass:.1e}\neps={soft:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
