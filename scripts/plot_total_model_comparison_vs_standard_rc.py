#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from horb_common import load_baryon_velocity_squared, has_einasto, einasto_label


def main():
    parser = argparse.ArgumentParser(
        description="Plot HORB and classical halo totals against a standard rotation curve."
    )
    parser.add_argument("--dm-compare-csv", required=True, type=Path)
    parser.add_argument("--rc-csv", required=True, type=Path)
    parser.add_argument("--baryons-csv", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--title", default="Total model comparison")
    parser.add_argument("--fit-report", type=Path, action="append", default=[])
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
    if has_einasto(dm):
        ax.plot(r, np.sqrt(dm["v_einasto_kms"]**2 + v_baryon2), linestyle="-.",
                linewidth=2.0, label=einasto_label(dm) + " + same baryons")
    ax.plot(r, v_baryon, linestyle=":", linewidth=1.8, label="literature baryons only")

    from horb_recorded import recorded_curve
    for report_path in args.fit_report:
        name,rr,v,report=recorded_curve(report_path,args.rc_csv,args.baryons_csv)
        ax.plot(rr,v,lw=2.5,label=f"{name} · RMS {report['rms_kms']:.2f} km/s")
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
