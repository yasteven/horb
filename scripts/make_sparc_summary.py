#!/usr/bin/env python3
"""Rebuild the README SPARC RMS and activity figures from checked snapshots.

--refit repeats the classical halo fits; it does not rerun the saved HORB search.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import differential_evolution, least_squares
from scipy.stats import spearmanr

from horb_highlights import G, halo_mass_shape, metrics
from horb_qi import velocity as qi_velocity, metadata as qi_metadata
from horb_vms_mond import fit_vms, mond_velocity, mond_metadata
from horb_sparc import load_data, baryonic_components, digest, save_json
from make_agn_catalog import build_catalog, concentration, load_sources, write_csv

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "data/sparc/density_comparison_20261126"
OUTPUT = ROOT / "data/sparc/model_comparison"
HALOS = ("nfw", "burkert", "piso", "einasto")
MODELS = ("squared", "psi", *HALOS, "vms", "mond", "qi", "baryons")
LABELS = {"squared": "HORB squared density", "psi": "HORB signed ψ",
          "nfw": "NFW", "burkert": "Burkert", "piso": "Pseudo-isothermal",
          "einasto": "Einasto", "vms": "VMS script profile", "mond": "MOND · simple",
          "qi": "Quantised inertia (QI)", "baryons": "Baryons alone"}
FIT_METHOD = {
    "version": 1, "objective": "sum(((v_model-v_observed)/published_error)^2)",
    "log10_mass_bounds_msun": [0., 14.], "log10_scale_bounds_kpc": [-3., 4.],
    "alpha_bounds": [.03, 2.], "reference_radius": "outermost measured radius",
    "optimizer": "differential_evolution followed by least_squares; explicit zero-halo candidate",
    "seed": 42, "maxiter": 600, "popsize": 12, "tol": 1e-10,
    "ml_disk": .5, "ml_bulge": .7,
}


def csv_hash(headers, rows):
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(headers)
    writer.writerows(rows)
    return hashlib.sha256(stream.getvalue().encode()).hexdigest()


def checked_inputs():
    catalog, curves, source = load_data()
    run = json.loads((SNAPSHOT / "run.json").read_text())
    if run["source"]["files"] != source["files"]:
        raise ValueError("HORB/SPARC source snapshots differ")
    if any(run["config"][k] != FIT_METHOD[k] for k in ("ml_disk", "ml_bulge")):
        raise ValueError("HORB/classical mass-to-light ratios differ")
    reports = json.loads((SNAPSHOT / "fit_reports.json").read_text())
    fits = {(r["galaxy"], r["density_model"]): r for r in reports}
    expected = {(g, m) for g in catalog for m in ("squared", "psi")}
    if len(fits) != len(reports) or set(fits) != expected:
        raise ValueError("expected exactly two HORB fits for each of 175 galaxies")
    for g, curve in curves.items():
        components = baryonic_components(curve, .5, .7)
        hashes = {
            "observations": csv_hash(["R_kpc", "Vobs_kms", "eV_kms"], curve[:, :3]),
            "baryons": csv_hash(["r_kpc", "v2_gas_kms2", "v2_disk_kms2", "v2_bulge_kms2", "v2_baryons_kms2"],
                               np.column_stack([curve[:, 0], components, components.sum(axis=1)])),
        }
        for model in ("squared", "psi"):
            fit = fits[g, model]
            if fit["n_observations"] != len(curve) or fit["config"]["criterion"] != "chi2":
                raise ValueError(f"{g}: wrong HORB observations or selection objective")
            for key, value in hashes.items():
                if fit["inputs"][key]["sha256"] != value:
                    raise ValueError(f"{g}: HORB and comparison {key} differ")
    return catalog, curves, fits, run


def halo_v2(model, radius, parameters, reference_radius):
    mass, scale = 10.**parameters[0], 10.**parameters[1]
    x = np.append(np.asarray(radius), reference_radius) / scale
    if model == "piso":
        # M(<r) is proportional to x-atan(x); avoid cancellation in the core.
        shape = np.where(x < .01, x**3*(1/3-x*x/5+x**4/7-x**6/9), x-np.arctan(x))
    else:
        alpha = 10.**parameters[2] if model == "einasto" else .18
        shape = halo_mass_shape(model, x, alpha)
    return G * mass / np.asarray(radius) * shape[:-1] / shape[-1]


def fit_halo(task):
    galaxy, model, curve = task
    radius, observed, errors = curve[:, :3].T
    baryons = baryonic_components(curve, .5, .7).sum(axis=1)
    bounds = [FIT_METHOD["log10_mass_bounds_msun"], FIT_METHOD["log10_scale_bounds_kpc"]]
    if model == "einasto":
        bounds.append(np.log10(FIT_METHOD["alpha_bounds"]).tolist())

    def residual(parameters):
        total = baryons + halo_v2(model, radius, parameters, radius[-1])
        if np.any(total < 0) or not np.isfinite(total).all():
            return np.full_like(observed, 1e50)
        return (np.sqrt(total)-observed)/errors

    global_fit = differential_evolution(lambda p: float(np.sum(residual(p)**2)), bounds,
        seed=FIT_METHOD["seed"], maxiter=FIT_METHOD["maxiter"], popsize=FIT_METHOD["popsize"],
        tol=FIT_METHOD["tol"], polish=False)
    local_fit = least_squares(residual, global_fit.x, bounds=np.array(bounds).T,
        jac="3-point", ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=3000)
    if not local_fit.success or np.any(np.abs(local_fit.fun) >= 1e49):
        raise ValueError(f"{galaxy}/{model}: local fit failed: {local_fit.message}")
    p = local_fit.x if np.sum(local_fit.fun**2) <= global_fit.fun else global_fit.x
    prediction = np.sqrt(baryons + halo_v2(model, radius, p, radius[-1]))
    zero = np.all(baryons >= 0) and np.sum(((np.sqrt(np.maximum(baryons, 0))-observed)/errors)**2) <= np.sum(residual(p)**2)
    if zero:
        prediction = np.sqrt(baryons)
    names = ["mass_within_reference_msun", "scale_kpc"] + (["alpha"] if model == "einasto" else [])
    return {"galaxy": galaxy, "model": model, "status": "ok", **metrics(prediction, observed, errors),
            "parameters": dict(zip(names, (10.**p).tolist())) if not zero else {"mass_within_reference_msun": 0.},
            "log10_parameters": p.tolist() if not zero else None, "zero_halo_selected": bool(zero),
            "reference_radius_kpc": float(radius[-1]), "parameter_count": len(p),
            "global_converged": bool(global_fit.success), "local_converged": bool(local_fit.success),
            "bounds_reached": [] if zero else [k for k, x, (lo, hi) in zip(names, p, bounds) if min(x-lo, hi-x) < 1e-4]}


def classical_fits(curves, refit=False, jobs=4):
    path = OUTPUT / "halo_fits.json"
    fingerprint = {"method": FIT_METHOD, "sparc_sources_sha256": digest(ROOT/"data/sparc/sources.json")}
    if refit:
        tasks = [(g, m, curves[g]) for g in sorted(curves) for m in HALOS]
        records = []
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            for i, fit in enumerate(pool.map(fit_halo, tasks), 1):
                records.append(fit)
                if i % 40 == 0 or i == len(tasks):
                    print(f"Classical halos: {i}/{len(tasks)}", flush=True)
        cache = {**fingerprint, "software": {"numpy": np.__version__, "scipy": scipy.__version__}, "fits": records}
        save_json(path, cache)
    else:
        if not path.exists():
            raise ValueError("missing classical fits; run with --refit")
        cache = json.loads(path.read_text())
        if any(cache.get(k) != v for k, v in fingerprint.items()):
            raise ValueError("classical fit cache is stale; run with --refit")
    records = cache["fits"]
    if len(records) != len(curves)*len(HALOS) or {(r["galaxy"], r["model"]) for r in records} != {(g,m) for g in curves for m in HALOS}:
        raise ValueError("classical fits are incomplete or duplicated")
    # Recompute every cached score from its physical parameters and original data.
    for record in records:
        curve = curves[record["galaxy"]]
        b = baryonic_components(curve, .5, .7).sum(axis=1)
        v2 = b if record["zero_halo_selected"] else b + halo_v2(record["model"], curve[:,0], record["log10_parameters"], curve[-1,0])
        if np.any(v2 < 0):
            raise ValueError("cached halo has outward total force")
        score = metrics(np.sqrt(v2), curve[:,1], curve[:,2])
        if not all(np.isclose(score[k], record[k], rtol=1e-9, atol=1e-9) for k in ("chi2", "rms_kms", "n")):
            raise ValueError("cached halo score does not match predictions")
    return records


def comparison_rows(curves, fits, halos):
    lookup = {(r["galaxy"], r["model"]): r for r in halos}
    result = []
    for galaxy, curve in sorted(curves.items()):
        b = baryonic_components(curve, .5, .7).sum(axis=1)
        for model in MODELS:
            if model in ("squared", "psi"):
                record = fits[galaxy, model]
                score = {k: record[k] for k in ("rms_kms", "chi2")}
                status, k = "ok", record["k"]
            elif model in HALOS:
                record = lookup[galaxy, model]
                score = {k: record[k] for k in ("rms_kms", "chi2")}
                status, k = "ok", record["parameter_count"]
            elif model == "vms":
                score = fit_vms(curve[:,0], curve[:,1], curve[:,2], b)
                status, k = "ok", 1
            elif np.any(b < 0):
                score, status, k = {"rms_kms": None, "chi2": None}, "undefined_outward_baryonic_force", 0
            else:
                prediction = (qi_velocity(curve[:,0], b) if model == "qi" else
                              mond_velocity(curve[:,0], b) if model == "mond" else np.sqrt(b))
                score = metrics(prediction, curve[:,1], curve[:,2])
                status, k = "ok", 0
            result.append({"galaxy": galaxy, "model": model, "status": status, "n_points": len(curve),
                           "rms_kms": score["rms_kms"], "chi2": score["chi2"], "fitted_coefficients_or_parameters": k})
    return result


def aggregate(rows):
    valid = {m: {r["galaxy"] for r in rows if r["model"] == m and r["status"] == "ok"} for m in MODELS}
    common = set.intersection(*valid.values())
    result = []
    for model in MODELS:
        subset = [r for r in rows if r["model"] == model and r["status"] == "ok"]
        rms = [r["rms_kms"] for r in subset]
        shared = [r["rms_kms"] for r in subset if r["galaxy"] in common]
        result.append({"model": model, "galaxies": len(rms), "mean_rms_kms": float(np.mean(rms)),
                       "median_rms_kms": float(np.median(rms)), "common_galaxies": len(common),
                       "common_mean_rms_kms": float(np.mean(shared))})
    return result


def activity_rows(fits):
    load_sources()  # Check hashes of radio tables, aliases and evidence annotations.
    rows = []
    for model in ("squared", "psi"):
        source_rows = {r["galaxy"]: r for r in build_catalog(density_model=model)[0] if r["catalog_coverage"] == "matched"}
        with (SNAPSHOT / f"agn_{model}.csv").open() as stream:
            matched = list(csv.DictReader(stream))
        if len(matched) != 51 or {r["galaxy"] for r in matched} != set(source_rows):
            raise ValueError("expected 51 distinct activity matches per density family")
        for row in matched:
            for key in ("activity_name", "optical_class", "radio_core_limit", "jet_evidence", "jet_evidence_source"):
                if row[key] != source_rows[row["galaxy"]][key]:
                    raise ValueError(f"activity snapshot differs from source catalog: {key}")
            if not np.isclose(float(row["log_core_Lnu_1p5_W_Hz_sparc_distance"]),
                              source_rows[row["galaxy"]]["log_core_Lnu_1p5_W_Hz_sparc_distance"], atol=1e-12, rtol=0):
                raise ValueError("activity luminosity differs from source catalog")
            fit = fits[row["galaxy"], model]
            measured = concentration(fit["components"])
            for key in ("largest_abs_mass_fraction", "effective_components"):
                if not np.isclose(float(row[key]), measured[key], rtol=1e-12):
                    raise ValueError("activity snapshot does not match selected HORB weights")
            optical = row["optical_class"]
            group = ("resolved_jets" if row["jet_evidence"] == "resolved_jet_lobes" else
                     "other_optical_AGN" if optical in ("Seyfert", "LINER") else
                     "HII_comparison" if optical == "HII" else "excluded_AGN_status_uncertain")
            rows.append({"galaxy": row["galaxy"], "density_model": model, "name": row["activity_name"],
                         "group": group, "optical_class": optical, "radio_core_limit": row["radio_core_limit"],
                         "log_Lnu_W_Hz": float(row["log_core_Lnu_1p5_W_Hz_sparc_distance"]),
                         "concentration": measured["largest_abs_mass_fraction"],
                         "effective_components": measured["effective_components"],
                         "coefficient_bound_reached": fit["coefficient_bound_reached"],
                         "rotation_rms_kms": fit["rms_kms"], "points": fit["n_observations"],
                         "sparc_quality": int(row["sparc_quality"]), "inclination_deg": float(row["inclination_deg"]),
                         "jet_evidence_source": row["jet_evidence_source"]})
    return rows


def activity_statistics(rows):
    result = {}
    for model in ("squared", "psi"):
        selected = [r for r in rows if r["density_model"] == model]
        agn = [r for r in selected if r["optical_class"] in ("Seyfert", "LINER")]
        detections = [r for r in agn if r["radio_core_limit"] == "detection"]
        rho, pvalue = spearmanr([r["log_Lnu_W_Hz"] for r in detections], [r["concentration"] for r in detections])
        groups = {}
        for group in ("resolved_jets", "other_optical_AGN", "HII_comparison"):
            members = [r for r in selected if r["group"] == group]
            groups[group] = {"n": len(members), "median_concentration": float(np.median([r["concentration"] for r in members])),
                             "coefficient_bound_count": sum(r["coefficient_bound_reached"] for r in members)}
        result[model] = {"optical_AGN": len(agn), "detections": len(detections), "upper_limits": len(agn)-len(detections),
                         "spearman_rho_detections_only": float(rho), "spearman_p_two_sided": float(pvalue), "groups": groups}
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refit", action="store_true", help="repeat all 700 classical halo fits")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    _, curves, fits, run = checked_inputs()
    halos = classical_fits(curves, args.refit, args.jobs)
    rows = comparison_rows(curves, fits, halos)
    summary = aggregate(rows)
    activity = activity_rows(fits)
    stats = activity_statistics(activity)
    write_csv(OUTPUT / "galaxy_scores.csv", rows)
    write_csv(OUTPUT / "mean_rms.csv", summary)
    write_csv(OUTPUT / "activity_concentrations.csv", activity)
    vms_fits = [{"galaxy":g, **fit_vms(c[:,0], c[:,1], c[:,2], baryonic_components(c,.5,.7).sum(axis=1))}
                for g,c in sorted(curves.items())]
    save_json(OUTPUT / "vms_fits.json", vms_fits)
    provenance_paths = [SNAPSHOT/"fit_reports.json", SNAPSHOT/"run.json", SNAPSHOT/"agn_squared.csv", SNAPSHOT/"agn_psi.csv",
                        ROOT/"data/sparc/sources.json", ROOT/"data/agn/sources.json", OUTPUT/"halo_fits.json", OUTPUT/"vms_fits.json"]
    save_json(OUTPUT / "summary.json", {"schema_version": 1, "galaxies_evaluated": len(curves),
        "rotation_points": sum(len(x) for x in curves.values()), "aggregation": "arithmetic mean of per-galaxy RMS; equal galaxy weights",
        "models": summary, "qi": qi_metadata(), "mond": mond_metadata(), "classical_method": FIT_METHOD, "horb_search": run["config"],
        "undefined_scores": [{k:r[k] for k in ("galaxy", "model", "status")} for r in rows if r["status"] != "ok"],
        "classical_bounds_reached": sum(bool(r["bounds_reached"]) for r in halos),
        "classical_global_nonconvergence": sum(not r["global_converged"] for r in halos),
        "activity_statistics": stats, "activity_exclusions": {"no_catalog_match": 124, "ALG_uncertain_AGN": 1},
        "input_sha256": {str(p.relative_to(ROOT)): digest(p) for p in provenance_paths}})
    from plot_sparc_summary import plot_rms, plot_activity
    plot_rms(summary)
    plot_activity(activity, stats)
    print(json.dumps({"mean_rms": summary, "activity": stats}, indent=2))


if __name__ == "__main__":
    main()
