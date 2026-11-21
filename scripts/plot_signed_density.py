#!/usr/bin/env python3
"""Show signed density without discarding negative regions."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv",type=Path)
    p.add_argument("-o","--out",required=True,type=Path)
    args=p.parse_args()
    data=np.atleast_1d(np.genfromtxt(args.csv,delimiter=",",names=True))
    x,z=np.unique(data["x_kpc"]),np.unique(data["z_kpc"])
    order=np.lexsort((data["x_kpc"],data["z_kpc"]))
    rho=data["rho_signed_msun_per_kpc3"][order].reshape(len(z),len(x))
    if not np.all(np.isfinite(rho)): raise ValueError("density must be finite")
    limit=max(float(np.max(np.abs(rho))),1e-30)
    fig,ax=plt.subplots(figsize=(9,7))
    im=ax.pcolormesh(x,z,rho,shading="auto",cmap="RdBu_r",norm=SymLogNorm(linthresh=limit*1e-4,vmin=-limit,vmax=limit))
    fig.colorbar(im,ax=ax,label="Signed density [M☉/kpc³] (symmetric log scale)")
    ax.set(xlabel="x [kpc]",ylabel="z [kpc]",title="Signed orbital-density sum: red positive, blue negative",aspect="equal")
    fig.tight_layout(); args.out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.out,dpi=180); plt.close(fig)
    print(f"wrote {args.out}")


if __name__=="__main__": main()
