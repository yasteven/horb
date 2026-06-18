#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

G_KPC = 4.3009e-6


def load_baryon_v2(path, r):
    rows = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
    if rows.shape == ():
        rows = np.array([rows], dtype=rows.dtype)

    v2 = np.zeros_like(r, dtype=float)

    for row in rows:
        kind = str(row["kind"])
        mass = float(row["mass_msun"])
        scale = float(row["scale_kpc"])

        if kind in ("disk", "exponential_disk", "exponential-disk"):
            x = r / scale
            menc = mass * (1.0 - np.exp(-x) * (1.0 + x))
        elif kind in ("bulge", "hernquist", "hernquist_bulge", "hernquist-bulge"):
            menc = mass * r**2 / (r + scale) ** 2
        else:
            raise ValueError(f"unknown baryon kind: {kind}")

        v = np.sqrt(np.where(r > 0, G_KPC * menc / r, 0.0))
        v2 += v * v

    return v2


def parse_coeffs(s):
    return [float(x) for x in str(s).split(";") if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-csv", required=True, type=Path)
    ap.add_argument("--rc-csv", required=True, type=Path)
    ap.add_argument("--baryons-csv", required=True, type=Path)
    ap.add_argument("--summary-out", required=True, type=Path)
    ap.add_argument("--bar-out", required=True, type=Path)
    ap.add_argument("--curve-out", required=True, type=Path)
    ap.add_argument("--metric", choices=["chi2_per_point", "rms_kms"], default="chi2_per_point")
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    scan = np.genfromtxt(args.scan_csv, delimiter=",", names=True, dtype=None, encoding=None)
    rc = np.genfromtxt(args.rc_csv, delimiter=",", names=True, dtype=None, encoding=None)

    if scan.shape == ():
        scan = np.array([scan], dtype=scan.dtype)

    order = np.argsort(scan[args.metric])
    ranked = scan[order]

    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.summary_out, "w") as f:
        f.write("rank,states,coeffs,a0_star_kpc,dm_mass_msun,softening_kpc,rms_kms,chi2_per_point\n")
        for i, row in enumerate(ranked[:100], start=1):
            f.write(
                f"{i},\"{row['states']}\",\"{row['coeffs']}\","
                f"{float(row['a0_star_kpc']):.6f},"
                f"{float(row['dm_mass_msun']):.6e},"
                f"{float(row['softening_kpc']):.6f},"
                f"{float(row['rms_kms']):.6f},"
                f"{float(row['chi2_per_point']):.6f}\n"
            )

    # Bar plot of top wavefunction candidates.
    top = ranked[: args.top]
    labels = []
    values = []

    for i, row in enumerate(top, start=1):
        coeffs = parse_coeffs(row["coeffs"])
        coeff_tag = ",".join(f"{c:.2f}" for c in coeffs)
        labels.append(
            f"#{i} a0={float(row['a0_star_kpc']):.2f} M={float(row['dm_mass_msun']):.1e}\n[{coeff_tag}]"
        )
        values.append(float(row[args.metric]))

    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    ax.bar(range(len(values)), values)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel(args.metric)
    ax.set_title(f"Top CUDA real-wavefunction candidates by {args.metric}")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    args.bar_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.bar_out, dpi=220)
    plt.close(fig)

    # Curve overlay: best wavefunction fit, using existing scan score only.
    # We cannot reconstruct the exact CUDA curve from scan CSV alone unless we rerun the candidate,
    # so this plot is a score/ranking artifact plus observed RC and baryons baseline.
    best = ranked[0]

    r = rc["R_kpc"]
    vb2 = load_baryon_v2(args.baryons_csv, r)

    fig, ax = plt.subplots(figsize=(11.0, 7.0))
    ax.errorbar(
        rc["R_kpc"],
        rc["Vobs_kms"],
        yerr=rc["eV_kms"],
        fmt="o",
        capsize=3,
        markersize=4,
        alpha=0.75,
        label="Sofue 2020 Milky Way RC",
    )
    ax.plot(r, np.sqrt(vb2), "--", linewidth=1.8, label="Literature baryons only")

    coeffs = parse_coeffs(best["coeffs"])
    coeff_label = ", ".join(f"{c:.3f}" for c in coeffs)

    text = (
        f"Best scanned real wavefunction\n"
        f"states = {best['states']}\n"
        f"coeffs = [{coeff_label}]\n"
        f"a0 = {float(best['a0_star_kpc']):.3f} kpc\n"
        f"M = {float(best['dm_mass_msun']):.3e} Msun\n"
        f"eps = {float(best['softening_kpc']):.3f} kpc\n"
        f"RMS = {float(best['rms_kms']):.3f} km/s\n"
        f"chi2/pt = {float(best['chi2_per_point']):.3f}"
    )

    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round", alpha=0.15),
    )

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v [km/s]")
    ax.set_title("Best CUDA real-wavefunction scan result")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    args.curve_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.curve_out, dpi=220)
    plt.close(fig)

    print(f"best chi2/pt={float(best['chi2_per_point']):.6f} rms={float(best['rms_kms']):.6f}")
    print(f"wrote {args.summary_out}")
    print(f"wrote {args.bar_out}")
    print(f"wrote {args.curve_out}")


if __name__ == "__main__":
    main()
