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

    # Rows were written z-major, then x-inner:
    # for z:
    #   for x:
    rho_grid = rho.reshape((nz, nx))

    return xs, zs, rho_grid


def plot_density(csv_path: Path, out_path: Path, title: str | None):
    xs, zs, rho = load_xz_csv(csv_path)

    # Avoid log10(0). The model should be positive except maybe exact nodes,
    # so use a tiny floor based on positive finite values.
    positive = rho[np.isfinite(rho) & (rho > 0.0)]
    if positive.size == 0:
        raise ValueError("No positive density values found.")

    floor = positive.min() * 1e-6
    log_rho = np.log10(np.maximum(rho, floor))

    fig, ax = plt.subplots(figsize=(8, 8))

    im = ax.imshow(
        log_rho,
        origin="lower",
        extent=[xs.min(), xs.max(), zs.min(), zs.max()],
        aspect="equal",
    )

    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("z [kpc]")

    if title is None:
        title = csv_path.stem

    ax.set_title(title)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10 density [M_sun / kpc^3]")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot HORB x-z density slice CSV as a log-density heatmap."
    )
    parser.add_argument("csv", type=Path, help="Input x-z density CSV.")
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=None,
        help="Output PNG path. Defaults to plots/<csv stem>.png",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Plot title. Defaults to CSV filename stem.",
    )

    args = parser.parse_args()

    csv_path = args.csv
    if args.out is None:
        out_path = Path("plots") / f"{csv_path.stem}.png"
    else:
        out_path = args.out

    out_path.parent.mkdir(parents=True, exist_ok=True)

    plot_density(csv_path, out_path, args.title)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
