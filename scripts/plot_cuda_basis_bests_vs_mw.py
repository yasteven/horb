#!/usr/bin/env python3

import argparse
from pathlib import Path
import subprocess
import os

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


def best_rows_by_state(scan, metric):
    states = sorted(set(str(s) for s in scan["state"]))

    bests = []
    for state in states:
        rows = scan[scan["state"] == state]
        idx = np.argmin(rows[metric])
        row = rows[idx]
        bests.append(row)

    bests.sort(key=lambda row: float(row[metric]))
    return bests


def run_curve_writer(row, lib_path):
    state = str(row["state"])
    a0 = float(row["a0_star_kpc"])
    mass = float(row["dm_mass_msun"])
    n_side = int(row["n_side"])
    extent = float(row["extent_kpc"])
    softening = float(row["softening_kpc"])

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{lib_path}:{env.get('LD_LIBRARY_PATH', '')}"

    cmd = [
        "cargo",
        "run",
        "-q",
        "-p",
        "cuda_kernels",
        "--bin",
        "write_horb_cuda_curve_basis",
        "--",
        state,
        str(a0),
        f"{mass:.12e}",
        str(n_side),
        str(extent),
        str(softening),
    ]

    result = subprocess.run(cmd, check=True, text=True, capture_output=True, env=env)
    csv_path = result.stdout.strip().splitlines()[-1].strip()

    return Path(csv_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-csv", required=True, type=Path)
    parser.add_argument("--rc-csv", required=True, type=Path)
    parser.add_argument("--baryons-csv", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--metric", choices=["chi2_per_point", "rms_kms"], default="chi2_per_point")
    parser.add_argument("--top", type=int, default=9)
    parser.add_argument("--lib-path", default="cuda_kernels/cuda")
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    scan = np.genfromtxt(args.scan_csv, delimiter=",", names=True, dtype=None, encoding=None)
    rc = np.genfromtxt(args.rc_csv, delimiter=",", names=True, dtype=None, encoding=None)

    bests = best_rows_by_state(scan, args.metric)[: args.top]

    fig, ax = plt.subplots(figsize=(11.5, 7.2))

    ax.errorbar(
        rc["R_kpc"],
        rc["Vobs_kms"],
        yerr=rc["eV_kms"],
        fmt="o",
        capsize=3,
        markersize=4,
        label="Sofue 2020 Milky Way RC",
        alpha=0.80,
    )

    print("Best rows:")
    print("state,a0_star_kpc,dm_mass_msun,n_side,extent_kpc,softening_kpc,rms_kms,chi2_per_point")

    for row in bests:
        state = str(row["state"])
        a0 = float(row["a0_star_kpc"])
        mass = float(row["dm_mass_msun"])
        n_side = int(row["n_side"])
        extent = float(row["extent_kpc"])
        soft = float(row["softening_kpc"])
        rms = float(row["rms_kms"])
        chi = float(row["chi2_per_point"])

        print(f"{state},{a0:.6f},{mass:.6e},{n_side},{extent:.6f},{soft:.6f},{rms:.6f},{chi:.6f}")

        curve_csv = run_curve_writer(row, args.lib_path)
        curve = np.genfromtxt(curve_csv, delimiter=",", names=True)

        r = curve["r_kpc"]
        vdm = curve["v_horb_cuda_kms"]
        vb2 = load_baryon_velocity_squared(args.baryons_csv, r)
        vtot = np.sqrt(vdm**2 + vb2)

        label = f"{state}: a0={a0:.2f}, M={mass:.1e}, eps={soft:.2f}, chi2/pt={chi:.2f}"
        ax.plot(r, vtot, linewidth=2.0, label=label)

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v [km/s]")
    title = args.title or f"Best CUDA disk-plane hydrogenic orbital fits vs Milky Way ({args.metric})"
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
