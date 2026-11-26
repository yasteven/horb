#!/usr/bin/env python3
"""Compare density definitions on identical data, including zero coefficients."""
import argparse
import json
from pathlib import Path
import numpy as np
from fit_signed_density import atomic_write, write_csv
from horb_highlights import metrics
from horb_lobes import MODEL_LABELS


def compare(prefixes, output, title, plot=True):
    reports = [json.loads(Path(str(p)+".json").read_text()) for p in prefixes]
    curves = [np.atleast_1d(np.genfromtxt(str(p)+"_curve.csv", delimiter=",", names=True)) for p in prefixes]
    observations = [np.atleast_1d(np.genfromtxt(str(p)+"_observations.csv", delimiter=",", names=True)) for p in prefixes]
    for report, obs in zip(reports[1:], observations[1:]):
        for k in ("observations", "baryons"):
            if report["inputs"][k]["sha256"] != reports[0]["inputs"][k]["sha256"]:
                raise ValueError("comparison inputs differ")
        for k in ("R_kpc", "Vobs_kms", "eV_kms", "used_in_fit"):
            if not np.array_equal(obs[k], observations[0][k]):
                raise ValueError("comparison observations or fit masks differ")
    obs = observations[0]
    used = obs["used_in_fit"].astype(bool)
    radius = curves[0]["r_kpc"]
    b = curves[0]["v2_baryons_kms2"]
    scores, speeds = {}, {}
    for report, curve, measured in zip(reports, curves, observations):
        model = report.get("density_model", "squared")
        if model in scores:
            raise ValueError("duplicate density model in comparison")
        if not np.array_equal(curve["r_kpc"], radius) or not np.allclose(curve["v2_baryons_kms2"], b):
            raise ValueError("comparison grids/baryons differ")
        score = metrics(measured["v_model_kms"][used], obs["Vobs_kms"][used], obs["eV_kms"][used])
        if not np.isclose(score["chi2"], report["chi2"], rtol=1e-8, atol=1e-8):
            raise ValueError("saved fit and comparison scores differ")
        scores[model] = {**score, "k":report["k"], "components":report["components"]}
        speeds[model] = curve["v_total_kms"]
    qb = np.interp(obs["R_kpc"][used], radius, b)
    scores["zero_orbitals"] = ({**metrics(np.sqrt(qb), obs["Vobs_kms"][used], obs["eV_kms"][used]),
                                "k":0, "components":[]} if np.all(qb >= 0) else
                               {"k":0, "components":[], "rms_kms":None, "chi2":None,
                                "status":"baryons alone have outward force at fitted radii"})
    speeds["zero_orbitals"] = np.sqrt(np.maximum(b, 0))
    speeds["zero_orbitals"][b < 0] = np.nan
    result = {"galaxy":title, "scores":scores, "inputs":reports[0]["inputs"],
              "selection":"same observations and fixed baryons; fitted models selected by their recorded criterion",
              "zero_orbitals":"all coefficients zero: baryons alone; Keplerian only outside finite enclosed baryonic mass"}
    output = str(output)
    atomic_write(output+".json", lambda f:json.dump(result, f, indent=2, allow_nan=False))
    write_csv(output+"_scores.csv", ["model","rms_kms","chi2","k"],
              [[m,s["rms_kms"],s["chi2"],s["k"]] for m,s in scores.items()])
    write_csv(output+"_curves.csv", ["r_kpc", *[f"v_{m}_kms" for m in speeds]],
              np.column_stack([radius,*speeds.values()]))
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(11,6))
        fig.subplots_adjust(right=.72)
        ax.errorbar(obs["R_kpc"], obs["Vobs_kms"], yerr=obs["eV_kms"], fmt="o", ms=3, color="black", label="Observations")
        colors = {"squared":"#138f92", "psi":"#bd6226", "zero_orbitals":"#858a95"}
        for i,(m,v) in enumerate(speeds.items()):
            ax.plot(radius, v, color=colors[m], lw=2.3, ls="--" if m=="zero_orbitals" else "-")
            rms=scores[m]["rms_kms"]
            name=MODEL_LABELS.get(m,"Zero orbitals / baryons")
            fig.text(.75,.79-i*.13,name,color=colors[m],size=12)
            fig.text(.75,.75-i*.13,"RMS undefined" if rms is None else f"RMS {rms:.3f} km/s",size=11)
        ax.set(title=title, xlabel="Radius [kpc]", ylabel="Circular speed [km/s]", xlim=(radius[0],radius[-1]))
        ax.grid(alpha=.2)
        fig.text(.75,.27,"Same fixed baryons\nScores on fitted data\nDifferent density bases",size=10)
        fig.savefig(output+".png",dpi=160);plt.close(fig)
    return result


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fit-prefix",type=Path,action="append",required=True)
    p.add_argument("--out-prefix",type=Path,required=True)
    p.add_argument("--title",default="HORB density comparison")
    p.add_argument("--no-plot",action="store_true")
    a=p.parse_args()
    compare(a.fit_prefix,a.out_prefix,a.title,not a.no_plot)
