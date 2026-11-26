"""Add provenance-checked fitted HORB families to fixed-parameter diagnostics."""
import json
from pathlib import Path
import numpy as np
from fit_signed_density import fingerprint
from horb_signed import baryons_function
from horb_swells.models import basis_matrix
from horb_lobes import MODEL_LABELS


def recorded_curve(report_path, observations_path, baryons_path):
    report=json.loads(Path(report_path).read_text())
    for key,path in (("observations",observations_path),("baryons",baryons_path)):
        if report["inputs"][key]["sha256"] != fingerprint(path)["sha256"]:
            raise ValueError(f"recorded HORB {key} differ; regenerate this fit before overlaying it")
    obs=np.atleast_1d(np.genfromtxt(observations_path,delimiter=",",names=True))
    lo,hi=report["r_fit_requested_kpc"]
    rr=obs["R_kpc"][(obs["R_kpc"]>=lo)&(obs["R_kpc"]<=hi)]
    r=np.unique(np.r_[np.linspace(lo,hi,800),rr])
    terms=report["components"]
    q=baryons_function(baryons_path)(r)+basis_matrix(r,terms)@np.array([t["weight"] for t in terms])
    if np.any(q<0): raise ValueError("recorded fit has negative total force")
    name=MODEL_LABELS[report.get("density_model","squared")]+" fitted"
    return name,r,np.sqrt(q),report
