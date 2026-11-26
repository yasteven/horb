#!/usr/bin/env python3
"""Reproduce both HORB definitions and all README baselines for AGC 114905."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from make_agc114905_data import prepare, ROOT, DATA


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--curve-bin",type=Path,required=True)
    p.add_argument("--refit",action="store_true")
    a=p.parse_args()
    prepare()
    obs=np.genfromtxt(DATA/"rotation_curve.csv",delimiter=",",names=True)
    lo,hi=obs["R_kpc"][[0,-1]]
    def run(command):
        subprocess.run([str(c) for c in command],check=True,cwd=ROOT)
    options=[]
    for model in ("squared","psi"):
        stem="agc114905"+("_psi" if model=="psi" else "")
        basis=ROOT/f"reports/{stem}_showcase_basis.csv"
        run([a.curve_bin,"signed_basis","--density-model",model,"--out",basis,"--scales","0.125,0.25,0.5,1,2",
             "--ref-mass","1e9","--r-min",lo,"--r-max",hi,"--r-step",(hi-lo)/256,"--rc-csv",DATA/"rotation_curve.csv"])
        report,weights=ROOT/f"data/highlights/{stem}_fit.json",ROOT/f"data/highlights/{stem}_weights.csv"
        prefix=ROOT/f"reports/{stem}_showcase"
        if a.refit or not report.exists():
            run([sys.executable,ROOT/"scripts/fit_signed_density.py","--basis-csv",basis,
                "--rc-csv",DATA/"rotation_curve.csv","--baryons-csv",DATA/"baryons.csv","--dataset-meta",DATA/"dataset.json",
                "--out-prefix",prefix,"--r-min",lo,"--r-max",hi,"--max-components","4","--search","exhaustive",
                "--budget","40000","--criterion","chi2","--ridge","0"])
            report,weights=Path(str(prefix)+".json"),Path(str(prefix)+"_weights.csv")
        density=ROOT/f"curves/{stem}_showcase_density_xz.csv"
        run([a.curve_bin,"signed_density","--weights",weights,"--out",density,"--extent","15","--bins","301"])
        flag="--psi-" if model=="psi" else "--"
        options.extend([flag+"fit-report",report,flag+"weights",weights,flag+"basis-csv",basis,flag+"density-csv",density])
    run([sys.executable,ROOT/"scripts/make_readme_highlights.py","--galaxy","agc114905",*options])
    record=json.loads((ROOT/"data/highlights/agc114905_comparison.json").read_text())
    print("Zero coefficients / baryons alone:",record["scores"]["baryons"])


if __name__ == "__main__":
    main()
