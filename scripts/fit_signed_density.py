#!/usr/bin/env python3
"""Fit a sparse signed linear combination of orbital densities to a rotation curve."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
import scipy
from horb_signed import load_basis, load_observations, baryons_function, Search


def atomic_write(path, write):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=path.name+".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="") as f:
            write(f)
        os.replace(temp, path)
    finally:
        if os.path.exists(temp): os.unlink(temp)


def write_csv(path, headers, rows):
    def write(f):
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)
    atomic_write(path, write)


def fingerprint(path):
    return {"path":str(Path(path).resolve()),"sha256":hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def dataset_metadata(path, rc_path, baryons_path):
    data = json.loads(Path(path).read_text())
    if data.get("schema_version") != 1 or not isinstance(data.get("galaxy"), str):
        raise ValueError("dataset metadata requires schema_version=1 and a galaxy name")
    for key, source in [("observations", rc_path), ("baryons", baryons_path)]:
        if data.get("inputs", {}).get(key, {}).get("sha256") != fingerprint(source)["sha256"]:
            raise ValueError(f"dataset metadata does not match {key}; regenerate the data or supply matching metadata")
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--basis-csv", required=True, type=Path)
    p.add_argument("--rc-csv", required=True, type=Path)
    p.add_argument("--baryons-csv", required=True, type=Path)
    p.add_argument("--dataset-meta", type=Path, help="dataset provenance JSON; its input hashes must match")
    p.add_argument("--out-prefix", required=True, type=Path)
    p.add_argument("--r-min", type=float, default=5)
    p.add_argument("--r-max", type=float, default=80)
    p.add_argument("--max-components", type=int, default=4)
    p.add_argument("--search", choices=["beam","anneal","exhaustive"], default="beam")
    p.add_argument("--beam-width", type=int, default=6)
    p.add_argument("--budget", type=int, default=2000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--max-abs-weight", type=float, default=20)
    p.add_argument("--ridge", type=float, default=1e-4)
    p.add_argument("--criterion", choices=["bic","chi2"], default="bic")
    p.add_argument("--weights", choices=["signed","nonnegative"], default="signed",
                   help="use nonnegative only as a control fit of the same basis")
    args = p.parse_args()
    provenance = dataset_metadata(args.dataset_meta, args.rc_csv, args.baryons_csv) if args.dataset_meta else None
    basis = load_basis(args.basis_csv)
    all_obs, fit_obs = load_observations(args.rc_csv,args.r_min,args.r_max,basis.radii)
    baryons = baryons_function(args.baryons_csv)
    rfit = fit_obs["R_kpc"]
    # Check positivity on every stored radius, intermediate subdivisions and
    # every fit observation. All interpolation is of signed force, never speed.
    rguard = np.unique(np.concatenate([basis.radii, rfit, *[
        basis.radii[:-1]*(1-t)+basis.radii[1:]*t for t in (0.25,0.5,0.75)]]))
    A, C = basis.at(rfit), basis.at(rguard)
    b, d = baryons(rfit), baryons(rguard)
    search = Search(A,b,fit_obs["Vobs_kms"],fit_obs["eV_kms"],C,d,
                    args.max_abs_weight,args.ridge,args.criterion,args.weights=="nonnegative")
    candidates = search.run(args.max_components,args.search,args.beam_width,args.budget,args.seed)
    best = candidates[0]
    terms = [{**basis.metadata[i],"weight":float(w),"signed_mass_msun":float(w*basis.metadata[i]["ref_mass_msun"])}
             for i,w in zip(best.subset,best.weights)]
    contributions = C[:,best.subset]*best.weights
    dark = np.sum(contributions,axis=1)
    total = d+dark
    if np.any(total<0): raise RuntimeError("internal error: invalid total force was accepted")
    prefix = str(args.out_prefix)
    weight_fields = ["basis_id","state","n","l","m","a0_star_kpc","ref_mass_msun","weight","signed_mass_msun"]
    write_csv(prefix+"_weights.csv",weight_fields,[[t[key] for key in weight_fields] for t in terms])
    write_csv(prefix+"_curve.csv",["r_kpc","v2_baryons_kms2","v2_dark_signed_kms2","v2_total_kms2","v_total_kms",
              *[f"v2_basis_{t['basis_id']}_kms2" for t in terms]],
              np.column_stack([rguard,d,dark,total,np.sqrt(total),contributions]))
    # Preserve every observation for the plot; mark the subset actually fitted.
    observed_rows = []
    for row in all_obs:
        radius = float(row["R_kpc"])
        supported = basis.radii[0] <= radius <= basis.radii[-1]
        pred = ""
        if supported:
            q = (baryons(np.array([radius])) + basis.at([radius])[:,best.subset]@best.weights).item()
            pred = np.sqrt(q) if q>=0 else ""
        used = args.r_min <= radius <= args.r_max
        observed_rows.append([radius,float(row["Vobs_kms"]),float(row["eV_kms"]),int(used),pred])
    write_csv(prefix+"_observations.csv",["R_kpc","Vobs_kms","eV_kms","used_in_fit","v_model_kms"],observed_rows)
    write_csv(prefix+"_candidates.csv",["rank","k","chi2","chi2_per_point","rms_kms","selection_score","basis_ids","weights","matrix_rank","condition_number"],
              [[rank,len(c.subset),c.chi2,c.chi2/len(rfit),c.rms,c.objective,
                ";".join(str(basis.metadata[i]["basis_id"]) for i in c.subset),
                ";".join(format(w,".17g") for w in c.weights),c.rank,c.condition]
               for rank,c in enumerate(candidates,1)])
    signed_masses = np.array([t["signed_mass_msun"] for t in terms])
    config = {key:(str(val) if isinstance(val,Path) else val) for key,val in vars(args).items()}
    report = {"schema_version":1,"model":"signed_linear_orbital_density","geometry":"azimuthal_mean",
              "config":config,"inputs":{"basis":fingerprint(args.basis_csv),"observations":fingerprint(args.rc_csv),"baryons":fingerprint(args.baryons_csv)},
              "software":{"numpy":np.__version__,"scipy":scipy.__version__},
              "basis_count":len(basis.metadata),"evaluated_subsets":len(search.cache),"failed_subsets":search.failed,
              "budget_reached":len(search.cache)>=args.budget,"n_observations":len(rfit),
              "r_fit_requested_kpc":[args.r_min,args.r_max],"r_guard_kpc":[float(rguard[0]),float(rguard[-1])],
              "guard_count":len(rguard),"k":len(terms),"chi2":best.chi2,"chi2_per_point":best.chi2/len(rfit),
              "rms_kms":best.rms,"selection_score":best.objective,"min_total_v2_kms2":best.min_total_v2,
              "matrix_rank":best.rank,"condition_number":best.condition if np.isfinite(best.condition) else None,
              "coefficient_bound_reached":bool(np.any(np.abs(best.weights)>=args.max_abs_weight*(1-1e-6))),
              "net_mass_msun":float(signed_masses.sum()),"positive_mass_msun":float(signed_masses[signed_masses>0].sum()),
              "negative_mass_msun":float(signed_masses[signed_masses<0].sum()),
              "components":terms,
              "selection_note":"chi2 + ridge*sum(w^2) + k*log(N) for bic; adaptive basis search makes this a heuristic complexity score, not a model probability",
              "domain_note":"positivity checked on the reported guard grid; no guarantee between samples or outside it; fit scores are in-sample"}
    if provenance is not None:
        report["dataset"] = provenance
        report["inputs"]["dataset_metadata"] = fingerprint(args.dataset_meta)
    # The JSON manifest is written last, after its associated tables succeed.
    atomic_write(prefix+".json",lambda f:json.dump(report,f,indent=2,allow_nan=False))
    print(f"best: k={len(terms)} RMS={best.rms:.6f} km/s chi2/N={best.chi2/len(rfit):.6f}")
    for term in terms:
        print(f"  {term['state']} a0={term['a0_star_kpc']:g} weight={term['weight']:+.8g} mass={term['signed_mass_msun']:+.8e} Msun")
    print(f"net mass={report['net_mass_msun']:.8e}; positive={report['positive_mass_msun']:.8e}; negative={report['negative_mass_msun']:.8e} Msun")
    print(f"evaluated {len(search.cache)} subsets; {search.failed} failed/infeasible; wrote {prefix}.json and CSV tables")


if __name__ == "__main__":
    main()
