#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def load_xz_csv(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)

    x = data["x_kpc"]
    z = data["z_kpc"]
    rho = data["rho_Msun_per_kpc3"]

    xs = np.unique(x)
    zs = np.unique(z)

    nx = len(xs)
    nz = len(zs)

    if nx * nz != len(rho):
        raise ValueError(
            f"CSV does not look like a complete grid: nx={nx}, nz={nz}, rows={len(rho)}"
        )

    rho_grid = rho.reshape((nz, nx))
    return xs, zs, rho_grid


def transformed_density(rho, mode: str, threshold: float):
    positive = rho[np.isfinite(rho) & (rho > 0.0)]
    if positive.size == 0:
        raise ValueError("No positive density values found.")

    rho_max = positive.max()

    if mode == "log":
        floor = positive.min() * 1e-6
        return np.log10(np.maximum(rho, floor)), "log10 density [M_sun / kpc^3]"

    if mode == "linear":
        return rho / rho_max, "normalized density"

    if mode == "threshold":
        norm = rho / rho_max
        masked = np.where(norm >= threshold, norm, np.nan)
        return masked, f"normalized density, threshold >= {threshold}"

    raise ValueError(f"unknown mode: {mode}")


def plot_density(csv_path: Path, out_path: Path, title: str | None, mode: str, threshold: float, contours: bool):
    xs, zs, rho = load_xz_csv(csv_path)
    image, cbar_label = transformed_density(rho, mode, threshold)

    fig, ax = plt.subplots(figsize=(8, 8))

    im = ax.imshow(
        image,
        origin="lower",
        extent=[xs.min(), xs.max(), zs.min(), zs.max()],
        aspect="equal",
    )

    if contours:
        norm = rho / np.nanmax(rho)
        levels = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75]
        ax.contour(xs, zs, norm, levels=levels, linewidths=0.8)

    ax.axhline(0.0, linewidth=0.8)
    ax.axvline(0.0, linewidth=0.8)

    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("z [kpc]")

    if title is None:
        title = f"{csv_path.stem} [{mode}]"

    ax.set_title(title)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot HORB x-z density slice CSV as a heatmap."
    )
    parser.add_argument("csv", type=Path, help="Input x-z density CSV.")
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=None,
        help="Output PNG path. Defaults to plots/<csv stem>_<mode>.png",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Plot title. Defaults to CSV filename stem.",
    )
    parser.add_argument(
        "--mode",
        choices=["log", "linear", "threshold"],
        default="log",
        help="Plot mode.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.15,
        help="Threshold for threshold mode, as fraction of max density.",
    )
    parser.add_argument(
        "--contours",
        action="store_true",
        help="Overlay normalized density contours.",
    )

    args = parser.parse_args()

    if args.out is None:
        args.out = Path("plots") / f"{args.csv.stem}_{args.mode}.png"

    args.out.parent.mkdir(parents=True, exist_ok=True)

    plot_density(
        csv_path=args.csv,
        out_path=args.out,
        title=args.title,
        mode=args.mode,
        threshold=args.threshold,
        contours=args.contours,
    )

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
