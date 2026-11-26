#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from horb_common import load_baryon_velocity_squared, has_einasto, einasto_label


def interp_to(r_target, r_source, v_source):
    return np.interp(r_target, r_source, v_source)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cuda-csv", required=True, type=Path)
    parser.add_argument("--dm-compare-csv", required=True, type=Path)
    parser.add_argument("--rc-csv", required=True, type=Path)
    parser.add_argument("--baryons-csv", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--title", default="Milky Way: CUDA HORB disk-plane force vs Sofue 2020")
    parser.add_argument("--fit-report", type=Path, action="append", default=[])
    args = parser.parse_args()

    cuda = np.genfromtxt(args.cuda_csv, delimiter=",", names=True)
    dm = np.genfromtxt(args.dm_compare_csv, delimiter=",", names=True)
    rc = np.genfromtxt(args.rc_csv, delimiter=",", names=True, dtype=None, encoding=None)

    common = (cuda["r_kpc"] >= dm["r_kpc"][0]) & (cuda["r_kpc"] <= dm["r_kpc"][-1])
    r = cuda["r_kpc"][common]
    v_cuda_dm = cuda["v_horb_cuda_kms"][common]

    v_baryon2 = load_baryon_velocity_squared(args.baryons_csv, r)
    v_baryon = np.sqrt(v_baryon2)

    v_cuda_total = np.sqrt(v_cuda_dm**2 + v_baryon2)

    v_sph_dm = interp_to(r, dm["r_kpc"], dm["v_horb_kms"])
    v_nfw_dm = interp_to(r, dm["r_kpc"], dm["v_nfw_kms"])
    v_piso_dm = interp_to(r, dm["r_kpc"], dm["v_piso_kms"])
    v_burkert_dm = interp_to(r, dm["r_kpc"], dm["v_burkert_kms"])

    v_sph_total = np.sqrt(v_sph_dm**2 + v_baryon2)
    v_nfw_total = np.sqrt(v_nfw_dm**2 + v_baryon2)
    v_piso_total = np.sqrt(v_piso_dm**2 + v_baryon2)
    v_burkert_total = np.sqrt(v_burkert_dm**2 + v_baryon2)

    fig, ax = plt.subplots(figsize=(10.5, 6.7))

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

    ax.plot(r, v_cuda_total, linewidth=3.0, label="HORB CUDA disk-plane + literature baryons")
    ax.plot(r, v_sph_total, linestyle="--", linewidth=2.2, label="HORB spherical approx + same baryons")
    ax.plot(r, v_nfw_total, linestyle="--", linewidth=1.9, label="NFW + same baryons")
    ax.plot(r, v_piso_total, linestyle="--", linewidth=1.9, label="pISO + same baryons")
    ax.plot(r, v_burkert_total, linestyle="--", linewidth=1.9, label="Burkert + same baryons")
    if has_einasto(dm):
        v_einasto = interp_to(r, dm["r_kpc"], dm["v_einasto_kms"])
        ax.plot(r, np.sqrt(v_einasto**2 + v_baryon2), linestyle="-.",
                linewidth=1.9, label=einasto_label(dm) + " + same baryons")
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
