#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Summarize HORB x-z density CSV.")
    parser.add_argument("csv", type=Path)
    args = parser.parse_args()

    data = np.genfromtxt(args.csv, delimiter=",", names=True)

    x = data["x_kpc"]
    z = data["z_kpc"]
    r = data["r_kpc"]
    theta = data["theta_rad"]
    rho = data["rho_Msun_per_kpc3"]

    imax = np.nanargmax(rho)

    positive = rho[np.isfinite(rho) & (rho > 0)]
    rho_max = positive.max()
    rho_min = positive.min()

    print(f"file: {args.csv}")
    print(f"rows: {len(rho)}")
    print(f"rho_min_positive: {rho_min:.6e}")
    print(f"rho_max: {rho_max:.6e}")
    print(f"log10_dynamic_range: {np.log10(rho_max / rho_min):.6f}")
    print()
    print("max density location:")
    print(f"  x_kpc: {x[imax]:.6f}")
    print(f"  z_kpc: {z[imax]:.6f}")
    print(f"  r_kpc: {r[imax]:.6f}")
    print(f"  theta_rad: {theta[imax]:.6f}")
    print()
    print("expected 3d_z2 spatial-density lobe peak:")
    print("  r_peak = 6 a0")
    print()
    print("expected 3d_z2 angular nodal cone:")
    print("  theta = acos(1/sqrt(3)) ≈ 0.955317 rad ≈ 54.7356 deg")
