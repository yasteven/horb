#!/usr/bin/env python3
"""Offline import of the attributed MP24 Figure 7 vector extraction.

Uses the exact published Figure 7 geometry, D=80.27 Mpc and i=35.86 deg,
as a common fixed scenario. Does not call this the fiducial measurement table.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from fit_signed_density import fingerprint, write_csv, atomic_write

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"data/agc114905"


def prepare(out=DATA):
    source = DATA/"mp24_figure7_vectors.json"
    meta = json.loads(source.read_text())
    axes = meta["axes"]
    x0, x2, y0, y10 = (axes[k] for k in ("x0_pdf","x2_pdf","y0_pdf","y10_pdf"))
    convert = lambda p: ((p[0]-x0)*2/(x2-x0), (y0-p[1])*10/(y0-y10))
    obs=[]
    for row in meta["observations"]:
        r, v = convert(row["point_pdf"])
        error = (row["error_bottom_pdf"]-row["error_top_pdf"])*5/(y0-y10)
        obs.append([r,v,float(np.hypot(error,.1)),error])
    obs=np.array(obs)
    if obs.shape != (5,4) or np.any(obs<=0) or np.any(np.diff(obs[:,0])<=0):
        raise ValueError("invalid digitized AGC 114905 observations")
    curves={k:np.array([convert(p) for p in points]) for k,points in meta["curve_polylines_pdf"].items()}
    r=np.unique(np.r_[np.linspace(obs[0,0],obs[-1,0],501),obs[:,0]])
    # Preserve the published total curve. Separate star/gas polylines have
    # different simplification nodes and need not sum exactly after interpolation.
    velocities={k:np.interp(r,c[:,0],c[:,1]) for k,c in curves.items()}
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/"rotation_curve.csv",["R_kpc","Vobs_kms","eV_kms","plotted_error_kms"],obs)
    write_csv(out/"baryons.csv",["r_kpc","v2_baryons_kms2"],np.column_stack([r,velocities["baryons"]**2]))
    write_csv(out/"digitized_baryonic_speeds.csv",["r_kpc","v_stars_kms","v_gas_kms","v_baryons_kms"],
              np.column_stack([r,velocities["stars"],velocities["gas"],velocities["baryons"]]))
    dataset={"schema_version":1,"galaxy":"AGC 114905","dataset_id":"agc114905_mp24_figure7_digitized",
        "source":meta,"vector_source":fingerprint(source),
        "rotation_curve":{"reference":meta["reference"]+" (digitized)","n":5},
        "baryons":{"plot_label":"Figure 7 stars + gas, fixed D=80.27 Mpc and i=35.86 deg"},
        "inputs":{"observations":fingerprint(out/"rotation_curve.csv"),"baryons":fingerprint(out/"baryons.csv")}}
    for entry in [dataset["vector_source"], *dataset["inputs"].values()]:
        path=Path(entry["path"])
        if path.is_relative_to(ROOT): entry["path"]=str(path.relative_to(ROOT))
    atomic_write(out/"dataset.json",lambda f:json.dump(dataset,f,indent=2,allow_nan=False))
    print("AGC 114905: five digitized 2024 points; Figure 7 geometry, not fiducial data; all-zero baseline retained.")
    return dataset


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir",type=Path,default=DATA)
    prepare(p.parse_args().out_dir)
