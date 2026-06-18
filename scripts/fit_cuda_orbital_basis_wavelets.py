#!/usr/bin/env python3

import argparse
import itertools
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


def nnls_projected_gradient(A, y, max_iter=20000, tol=1e-10):
    # Minimize ||A w - y||^2 with w >= 0.
    # Small/simple dependency-free NNLS.
    m, n = A.shape
    w = np.zeros(n)

    gram = A.T @ A
    rhs = A.T @ y

    L = np.linalg.norm(gram, 2)
    if L <= 0:
        return w

    step = 1.0 / L

    prev = np.inf

    for _ in range(max_iter):
        grad = gram @ w - rhs
        w = np.maximum(0.0, w - step * grad)

        obj = np.mean((A @ w - y) ** 2)
        if abs(prev - obj) < tol * max(1.0, obj):
            break
        prev = obj

    return w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-csv", required=True, type=Path)
    parser.add_argument("--rc-csv", required=True, type=Path)
    parser.add_argument("--baryons-csv", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    parser.add_argument("--summary-out", required=True, type=Path)
    parser.add_argument("--r-min", type=float, default=5.0)
    parser.add_argument("--r-max", type=float, default=25.0)
    parser.add_argument("--max-components", type=int, default=4)
    parser.add_argument("--top-basis", type=int, default=18)
    parser.add_argument("--metric", choices=["rms", "chi2"], default="chi2")
    args = parser.parse_args()

    basis = np.genfromtxt(args.basis_csv, delimiter=",", names=True, dtype=None, encoding=None)
    rc = np.genfromtxt(args.rc_csv, delimiter=",", names=True, dtype=None, encoding=None)

    mask = (rc["R_kpc"] >= args.r_min) & (rc["R_kpc"] <= args.r_max)
    r_fit = rc["R_kpc"][mask]
    vobs = rc["Vobs_kms"][mask]
    ev = rc["eV_kms"][mask]

    vb2 = load_baryon_v2(args.baryons_csv, r_fit)
    y = vobs**2 - vb2

    # Do not let negative target dark contribution blow up the fit.
    y = np.maximum(y, 0.0)

    basis_ids = sorted(set(int(x) for x in basis["basis_id"]))

    curves = []
    labels = []

    for bid in basis_ids:
        rows = basis[basis["basis_id"] == bid]
        state = str(rows["state"][0])
        a0 = float(rows["a0_star_kpc"][0])
        ref_mass = float(rows["ref_mass_msun"][0])

        r_curve = rows["r_kpc"].astype(float)
        v_curve = rows["v_cuda_kms"].astype(float)

        v_interp = np.interp(r_fit, r_curve, v_curve)
        dark_v2_per_ref_mass = v_interp**2

        # Score each basis alone with optimal nonnegative amplitude.
        denom = np.dot(dark_v2_per_ref_mass, dark_v2_per_ref_mass)
        amp = 0.0 if denom <= 0 else max(0.0, np.dot(dark_v2_per_ref_mass, y) / denom)
        pred = amp * dark_v2_per_ref_mass

        rms = np.sqrt(np.mean((np.sqrt(pred + vb2) - vobs) ** 2))
        chi2 = np.mean(((np.sqrt(pred + vb2) - vobs) / ev) ** 2)

        curves.append({
            "basis_id": bid,
            "state": state,
            "a0": a0,
            "ref_mass": ref_mass,
            "curve": dark_v2_per_ref_mass,
            "single_amp": amp,
            "single_rms": rms,
            "single_chi2": chi2,
        })

        labels.append(f"{state} a0={a0:.2f}")

    rank_key = "single_chi2" if args.metric == "chi2" else "single_rms"
    curves = sorted(curves, key=lambda x: x[rank_key])[: args.top_basis]

    best = None
    results = []

    for k in range(1, args.max_components + 1):
        for combo in itertools.combinations(range(len(curves)), k):
            A = np.column_stack([curves[i]["curve"] for i in combo])
            w = nnls_projected_gradient(A, y)
            pred_dark_v2 = A @ w
            vtot = np.sqrt(pred_dark_v2 + vb2)

            rms = np.sqrt(np.mean((vtot - vobs) ** 2))
            chi2 = np.mean(((vtot - vobs) / ev) ** 2)

            result = {
                "k": k,
                "combo": combo,
                "weights": w,
                "rms": rms,
                "chi2": chi2,
            }

            results.append(result)

            score = chi2 if args.metric == "chi2" else rms
            if best is None or score < best["score"]:
                best = {**result, "score": score}

    args.summary_out.parent.mkdir(parents=True, exist_ok=True)

    with open(args.summary_out, "w") as f:
        f.write("rank,k,rms_kms,chi2_per_point,total_dm_mass_msun,components\n")

        results_sorted = sorted(results, key=lambda x: x["chi2" if args.metric == "chi2" else "rms"])

        for rank, res in enumerate(results_sorted[:50], start=1):
            parts = []
            total_mass = 0.0

            for local_idx, weight in zip(res["combo"], res["weights"]):
                c = curves[local_idx]
                mass = weight * c["ref_mass"]
                total_mass += mass
                parts.append(f"{c['state']}@a0={c['a0']:.3f}:w={weight:.6g}:M={mass:.6e}")

            f.write(
                f"{rank},{res['k']},{res['rms']:.6f},{res['chi2']:.6f},{total_mass:.6e},"
                + "\""
                + "; ".join(parts)
                + "\"\n"
            )

    print("best:")
    print(f"k={best['k']} rms={best['rms']:.6f} chi2/pt={best['chi2']:.6f}")

    for local_idx, weight in zip(best["combo"], best["weights"]):
        c = curves[local_idx]
        print(
            f"  {c['state']} a0={c['a0']:.3f} weight={weight:.6g} "
            f"mass={weight * c['ref_mass']:.6e}"
        )

    # Plot best curve.
    r_plot = np.linspace(0.5, 80.0, 160)
    vb2_plot = load_baryon_v2(args.baryons_csv, r_plot)

    pred_dark_v2_plot = np.zeros_like(r_plot)

    for local_idx, weight in zip(best["combo"], best["weights"]):
        c = curves[local_idx]
        rows = basis[basis["basis_id"] == c["basis_id"]]
        r_curve = rows["r_kpc"].astype(float)
        v_curve = rows["v_cuda_kms"].astype(float)
        pred_dark_v2_plot += weight * np.interp(r_plot, r_curve, v_curve) ** 2

    vtot_plot = np.sqrt(pred_dark_v2_plot + vb2_plot)

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

    ax.plot(
        r_plot,
        vtot_plot,
        linewidth=2.8,
        label=f"Best sparse orbital-basis fit: k={best['k']}, RMS={best['rms']:.2f}, chi2/pt={best['chi2']:.2f}",
    )

    ax.plot(
        r_plot,
        np.sqrt(vb2_plot),
        linestyle="--",
        linewidth=1.5,
        label="Literature baryons only",
    )

    ax.set_xlabel("R [kpc]")
    ax.set_ylabel("v [km/s]")
    ax.set_title("Sparse CUDA hydrogenic orbital-basis fit to Milky Way rotation curve")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220)
    plt.close(fig)

    print(f"wrote {args.out}")
    print(f"wrote {args.summary_out}")


if __name__ == "__main__":
    main()
