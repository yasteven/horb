#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

G_KPC = 4.3009e-6


def baryon_velocity_component(kind, mass, scale, r):
    r = np.asarray(r)

    if kind in ("disk", "exponential_disk", "exponential-disk"):
        x = r / scale
        menc = mass * (1.0 - np.exp(-x) * (1.0 + x))
    elif kind in ("bulge", "hernquist", "hernquist_bulge", "hernquist-bulge"):
        menc = mass * r**2 / (r + scale) ** 2
    else:
        raise ValueError(f"unknown baryon component kind: {kind}")

    return np.sqrt(np.where(r > 0, G_KPC * menc / r, 0.0))


def load_baryon_velocity_squared(path, r):
    rows = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)

    if rows.shape == ():
        rows = np.array([rows], dtype=rows.dtype)

    v2 = np.zeros_like(r, dtype=float)

    for row in rows:
        kind = str(row["kind"])
        mass = float(row["mass_msun"])
        scale = float(row["scale_kpc"])
        v = baryon_velocity_component(kind, mass, scale, r)
        v2 += v**2

    return v2


def main():
    parser = argparse.ArgumentParser(
        description="Plot HORB and classical halo totals against a standard rotation curve."
    )
    parser.add_argument("--dm-compare-csv", required=True, type=Path)
    parser.add_argument("--rc-csv", required=True, type=Path)
    parser.add_argument("--baryons-csv", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--title", default="Total model comparison")
    args = parser.parse_args()

    dm = np.genfromtxt(args.dm_compare_csv, delimiter=",", names=True)
    rc = np.genfromtxt(args.rc_csv, delimiter=",", names=True, dtype=None, encoding=None)

    r = dm["r_kpc"]
    v_baryon2 = load_baryon_velocity_squared(args.baryons_csv, r)
    v_baryon = np.sqrt(v_baryon2)

    v_horb_total = np.sqrt(dm["v_horb_kms"]**2 + v_baryon2)
    v_piso_total = np.sqrt(dm["v_piso_kms"]**2 + v_baryon2)
    v_nfw_total = np.sqrt(dm["v_nfw_kms"]**2 + v_baryon2)
    v_burkert_total = np.sqrt(dm["v_burkert_kms"]**2 + v_baryon2)

    fig, ax = plt.subplots(figsize=(10, 6.5))

    ax.errorbar(
        rc["R_kpc"],
        rc["Vobs_kms"],
        yerr=rc["eV_kms"],
        fmt="o",
        capsize=3,
        markersize=4,
        label="Sofue 2020 Milky Way RC",
        alpha=0.85,
    )

    ax.plot(r, v_horb_total, linewidth=2.8, label="HORB 3d_z2 + literature baryons")
    ax.plot(r, v_nfw_total, linestyle="--", linewidth=2.0, label="NFW + same baryons")
    ax.plot(r, v_piso_total, linestyle="--", linewidth=2.0, label="pISO + same baryons")
    ax.plot(r, v_burkert_total, linestyle="--", linewidth=2.0, label="Burkert + same baryons")
    ax.plot(r, v_baryon, linestyle=":", linewidth=1.8, label="literature baryons only")

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
