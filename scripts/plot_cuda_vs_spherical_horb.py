#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cuda-csv", required=True, type=Path)
    parser.add_argument("--spherical-csv", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    args = parser.parse_args()

    cuda = np.genfromtxt(args.cuda_csv, delimiter=",", names=True)
    sph = np.genfromtxt(args.spherical_csv, delimiter=",", names=True)

    fig, ax = plt.subplots(figsize=(9.5, 6))

    ax.plot(
        cuda["r_kpc"],
        cuda["v_horb_cuda_kms"],
        linewidth=2.8,
        label="HORB actual disk-plane CUDA",
    )

    # compare_dm_fixed_baselines output has v_horb_kms
    ax.plot(
        sph["r_kpc"],
        sph["v_horb_kms"],
        linestyle="--",
        linewidth=2.2,
        label="HORB spherical enclosed-mass approximation",
    )

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v_DM [km/s]")
    ax.set_title("HORB 3d_z2: actual disk-plane force vs spherical approximation")
    ax.grid(True, alpha=0.3)
    ax.legend()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=220)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
