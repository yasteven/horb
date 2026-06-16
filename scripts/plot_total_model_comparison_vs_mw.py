#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Plot HORB+baryons and classical halo+baryons against Milky Way rotation curve."
    )
    parser.add_argument("--total-csv", required=True, type=Path)
    parser.add_argument("--dm-compare-csv", required=True, type=Path)
    parser.add_argument("--mw-csv", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--title", default="Total rotation-curve comparison vs Milky Way")
    args = parser.parse_args()

    total = np.genfromtxt(args.total_csv, delimiter=",", names=True)
    dm = np.genfromtxt(args.dm_compare_csv, delimiter=",", names=True)
    mw = np.genfromtxt(args.mw_csv, delimiter=",", names=True, dtype=None, encoding=None)

    r = total["r_kpc"]

    # Baryons from the best-candidate total CSV.
    v_disk = total["v_disk_kms"]
    v_bulge = total["v_bulge_kms"]
    v_baryon2 = v_disk**2 + v_bulge**2

    # Interpolate DM comparison onto the total-curve radius grid if needed.
    v_horb_dm = np.interp(r, dm["r_kpc"], dm["v_horb_kms"])
    v_piso_dm = np.interp(r, dm["r_kpc"], dm["v_piso_kms"])
    v_nfw_dm = np.interp(r, dm["r_kpc"], dm["v_nfw_kms"])
    v_burkert_dm = np.interp(r, dm["r_kpc"], dm["v_burkert_kms"])

    v_horb_total = np.sqrt(v_horb_dm**2 + v_baryon2)
    v_piso_total = np.sqrt(v_piso_dm**2 + v_baryon2)
    v_nfw_total = np.sqrt(v_nfw_dm**2 + v_baryon2)
    v_burkert_total = np.sqrt(v_burkert_dm**2 + v_baryon2)
    v_baryon_total = np.sqrt(v_baryon2)

    fig, ax = plt.subplots(figsize=(10, 6.5))

    ax.errorbar(
        mw["R_kpc"],
        mw["v_kms"],
        yerr=mw["v_err_kms"],
        fmt="o",
        capsize=3,
        markersize=4,
        label="Milky Way target",
        alpha=0.85,
    )

    ax.plot(r, v_horb_total, linewidth=2.8, label="HORB 3d_z2 + baryons")
    ax.plot(r, v_nfw_total, linestyle="--", linewidth=2.0, label="NFW + baryons")
    ax.plot(r, v_piso_total, linestyle="--", linewidth=2.0, label="pISO + baryons")
    ax.plot(r, v_burkert_total, linestyle="--", linewidth=2.0, label="Burkert + baryons")

    ax.plot(r, v_baryon_total, linestyle=":", linewidth=1.8, label="baryons only")

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v [km/s]")
    ax.set_title(args.title)
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
