#!/usr/bin/env python3
"""Plot the fitted rotation curve and the signed force contributions."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fit-prefix",required=True)
    p.add_argument("-o","--out",required=True,type=Path)
    args=p.parse_args()
    report=json.loads(Path(args.fit_prefix+".json").read_text())
    curve=np.atleast_1d(np.genfromtxt(args.fit_prefix+"_curve.csv",delimiter=",",names=True))
    obs=np.atleast_1d(np.genfromtxt(args.fit_prefix+"_observations.csv",delimiter=",",names=True))
    fig,axes=plt.subplots(2,1,figsize=(11,9),sharex=True)
    r=curve["r_kpc"]
    used=obs["used_in_fit"]==1
    for mask,label,color in [(used,"Fitted observations","black"),(~used,"Outside fit range","0.65")]:
        if np.any(mask): axes[0].errorbar(obs["R_kpc"][mask],obs["Vobs_kms"][mask],yerr=obs["eV_kms"][mask],fmt="o",ms=3,color=color,label=label)
    axes[0].plot(r,curve["v_total_kms"],label=("Signed ψ lobes + baryons" if report.get("density_model") == "psi" else "Squared density + baryons"),lw=2.5)
    qb=curve["v2_baryons_kms2"]
    axes[0].plot(r,np.sqrt(np.where(qb>=0,qb,np.nan)),"--",label="Fixed baryonic model")
    axes[0].set_ylabel("Circular speed [km/s]")
    dataset=report.get("dataset",{})
    galaxy=dataset.get("galaxy", "Signed orbital-density fit")
    axes[0].set_title(f"{galaxy}: k={report['k']}, RMS={report['rms_kms']:.2f} km/s, chi²/N={report['chi2_per_point']:.2f}")
    for term in report["components"]:
        axes[1].plot(r,curve[f"v2_basis_{term['basis_id']}_kms2"],label=f"{term['state']} a₀={term['a0_star_kpc']:g}, w={term['weight']:+.4g}")
    axes[1].plot(r,curve["v2_dark_signed_kms2"],color="black",lw=2,label="Sum of dark contributions")
    axes[1].axhline(0,color="0.5",lw=0.8)
    axes[1].set_ylabel("Signed R × inward acceleration [(km/s)²]")
    axes[1].set_xlabel("R [kpc]")
    for ax in axes:
        ax.axvspan(*report["r_fit_requested_kpc"],alpha=0.06,color="green")
        ax.grid(alpha=0.25); ax.legend(fontsize=8)
    if dataset:
        low,high=report["r_fit_requested_kpc"]
        margin=(high-low)*0.07
        axes[1].set_xlim(max(r[0],low-margin),min(r[-1],high+margin))
        source=dataset.get("rotation_curve",{}).get("reference","Supplied rotation curve")
        baryon_label=dataset.get("baryons",{}).get("plot_label","fixed baryonic model")
        domain_note = ("Published radii only; fixed baryons; in-sample fit."
                       if dataset.get("dataset_id") == "sparc_lelli2016" else
                       "Fit evaluated only in the shaded range; curves outside it are extrapolations.")
        note=source+" · "+baryon_label+"\n"+domain_note
        fig.text(0.5,0.015,note,ha="center",fontsize=8,color="0.35")
    fig.tight_layout(rect=(0,0.06 if dataset else 0,1,1))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.out,dpi=180); plt.close(fig)
    print(f"wrote {args.out}")


if __name__=="__main__": main()
