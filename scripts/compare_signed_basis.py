#!/usr/bin/env python3
"""Compare signed force libraries to check quadrature/CUDA convergence."""
import argparse
import csv
import sys
import numpy as np
from horb_signed import load_basis


def compare(reference,candidate):
    if len(reference.metadata)!=len(candidate.metadata):
        raise ValueError("different basis counts")
    r=reference.radii
    actual=candidate.at(r)
    rows=[]
    for i,(a,b) in enumerate(zip(reference.metadata,candidate.metadata)):
        for key in ("n","l","m","a0_star_kpc","ref_mass_msun","geometry","normalization","softening_kpc"):
            if a[key]!=b[key]: raise ValueError(f"basis {i}: different {key}")
        expected=reference.matrix[:,i]
        diff=actual[:,i]-expected
        denom=np.linalg.norm(expected)
        if denom==0: raise ValueError("zero reference force norm")
        rows.append([a["basis_id"],a["state"],a["a0_star_kpc"],float(np.linalg.norm(diff)/denom),float(np.max(np.abs(diff)))])
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",required=True)
    p.add_argument("--candidate",required=True)
    p.add_argument("--tolerance",type=float,default=0.01,help="maximum relative L2 error per basis (default 1%%)")
    args=p.parse_args()
    if not np.isfinite(args.tolerance) or args.tolerance<=0: p.error("tolerance must be finite and positive")
    rows=compare(load_basis(args.reference),load_basis(args.candidate))
    writer=csv.writer(sys.stdout,lineterminator="\n")
    writer.writerow(["basis_id","state","a0_star_kpc","relative_l2","max_abs_v2_kms2"])
    writer.writerows(rows)
    if any(row[3]>args.tolerance for row in rows):
        print("force difference exceeds tolerance; increase quadrature resolution before interpreting fitted weights",file=sys.stderr)
        return 1
    return 0


if __name__=="__main__": sys.exit(main())
