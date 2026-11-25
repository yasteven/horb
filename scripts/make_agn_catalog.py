#!/usr/bin/env python3
"""Join measured nuclear activity to SPARC and its existing HORB fits, offline."""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from horb_sparc import ROOT, baryonic_components, load_data
from fit_signed_density import atomic_write

DATA = ROOT / "data/agn"


def canonical(name):
    name = re.sub(r"\s", "", name.upper())
    match = re.fullmatch(r"(NGC|IC|UGC|UGCA|DDO|PGC)0*(\d+)", name)
    return match[1] + str(int(match[2])) if match else name


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"empty output: {path}")
    def write(stream):
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    atomic_write(Path(path), write)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_sources(directory=DATA):
    directory = Path(directory)
    manifest = json.loads((directory / "sources.json").read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported AGN catalog schema")
    for name, metadata in manifest["files"].items():
        if Path(name).name != name or sha256(directory / name) != metadata["sha256"]:
            raise ValueError(f"AGN source checksum mismatch: {name}")
    return manifest


def index_unique(rows, key):
    result = {}
    for row in rows:
        name = canonical(row[key])
        if name in result:
            raise ValueError(f"duplicate catalog name: {name}")
        result[name] = row
    return result


def concentration(weights):
    """Fractions of absolute signed mass, not probabilities or enclosed mass."""
    values = [float(row["signed_mass_msun"]) for row in weights]
    if not values or not all(math.isfinite(x) for x in values):
        raise ValueError("weights must contain finite signed masses")
    total = math.fsum(abs(x) for x in values)
    if total == 0:
        return {key: "" for key in (
            "largest_abs_mass_fraction", "effective_components", "weight_entropy_nats",
            "negative_abs_mass_fraction", "mass_cancellation_fraction")}
    fractions = [abs(x) / total for x in values]
    return {
        "largest_abs_mass_fraction": max(fractions),
        "effective_components": 1 / math.fsum(x*x for x in fractions),
        "weight_entropy_nats": -math.fsum(x*math.log(x) for x in fractions if x),
        "negative_abs_mass_fraction": math.fsum(abs(x) for x in values if x < 0) / total,
        "mass_cancellation_fraction": max(0.0, 1 - abs(math.fsum(values)) / total),
    }


def validate_weights(summary, weights):
    if len(weights) != int(summary["k"]):
        raise ValueError(f"{summary['galaxy']}: fit component count disagrees")
    seen = set()
    for row in weights:
        key = (row["state"], float(row["a0_star_kpc"]))
        if key in seen:
            raise ValueError("duplicate orbital component")
        seen.add(key)
        w, ref, mass = (float(row[k]) for k in ("weight", "ref_mass_msun", "signed_mass_msun"))
        if not all(map(math.isfinite, (w, ref, mass, key[1]))) or ref <= 0 or key[1] <= 0:
            raise ValueError("invalid orbital weight/scale")
        if not math.isclose(w*ref, mass, rel_tol=1e-10, abs_tol=1e-5):
            raise ValueError("signed mass disagrees with coefficient normalization")
    for key, values in (
        ("net_mass_msun", [float(r["signed_mass_msun"]) for r in weights]),
        ("positive_mass_msun", [float(r["signed_mass_msun"]) for r in weights if float(r["weight"]) > 0]),
        ("negative_mass_msun", [float(r["signed_mass_msun"]) for r in weights if float(r["weight"]) < 0]),
    ):
        if not math.isclose(math.fsum(values), float(summary[key]), rel_tol=1e-9, abs_tol=1e-3):
            raise ValueError(f"{summary['galaxy']}: summary and weights disagree on {key}")


def load_fits(fits_dir=None):
    weights = defaultdict(list)
    if fits_dir is None:
        summary_path = ROOT / "data/sparc/reference_summary.csv"
        run_path = ROOT / "data/sparc/reference_run.json"
        for row in read_csv(DATA / "reference_fit_weights.csv"):
            weights[row["galaxy"]].append(row)
        inputs = {str(p.relative_to(ROOT)): sha256(p) for p in
                  (summary_path, run_path, DATA / "reference_fit_weights.csv")}
        manifest = json.loads((DATA / "sources.json").read_text())
        if (sha256(summary_path) != manifest["fit_summary_sha256"]
                or sha256(run_path) != manifest["fit_run_sha256"]):
            raise ValueError("bundled SPARC reference fit snapshot has changed")
    else:
        fits_dir = Path(fits_dir)
        summary_path, run_path = fits_dir / "summary.csv", fits_dir / "run.json"
        inputs = {str(p): sha256(p) for p in (summary_path, run_path)}
    summaries = index_unique(read_csv(summary_path), "galaxy")
    for summary in summaries.values():
        name = summary["galaxy"]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
            raise ValueError("invalid fit galaxy name")
        if summary["status"] != "complete":
            continue
        if fits_dir is not None:
            path = fits_dir / name / "fit_weights.csv"
            weights[name] = read_csv(path)
            inputs[str(path)] = sha256(path)
        validate_weights(summary, weights[name])
    return summaries, weights, json.loads(run_path.read_text()), inputs


def build_catalog(fits_dir=None):
    manifest = load_sources()
    catalog, curves, sparc_source = load_data()
    fits, weights, run, fit_inputs = load_fits(fits_dir)
    if run["source"]["files"] != sparc_source["files"]:
        raise ValueError("fit run uses a different SPARC source snapshot")
    radio = index_unique(read_csv(DATA / "lemmings_radio_2021.csv"), "canonical_name")
    xray = index_unique(read_csv(DATA / "lemmings_xray_2022.csv"), "canonical_name")
    aliases = index_unique(read_csv(DATA / "simbad_aliases.csv"), "id")
    evidence = {row["sparc_name"]: row for row in read_csv(DATA / "jet_evidence.csv")}
    rows, joined_weights, points = [], [], []
    ml_disk, ml_bulge = float(run["config"]["ml_disk"]), float(run["config"]["ml_bulge"])
    metric_names = list(concentration([{"signed_mass_msun": 1}]))
    for name, meta in sorted(catalog.items()):
        key = canonical(name)
        keys = {key} | {canonical(n) for n in aliases.get(key, {}).get("ids", "").split("|")}
        matches = keys & radio.keys()
        if len(matches) > 1:
            raise ValueError(f"ambiguous galaxy alias: {name}")
        match = next(iter(matches), None)
        r, x = radio.get(match, {}), xray.get(match, {})
        f = fits.get(key, {})
        curve = curves[name]
        activity_class = r.get("optical_class", "unclassified")
        core_detection = r.get("log_core_nuLnu_erg_s_limit") == "detection"
        jet = evidence.get(name, {})
        jet_class = jet.get("jet_evidence", "unassessed")
        if match and not jet:
            jet_class = ("compact_radio_AGN_candidate" if core_detection and activity_class in ("Seyfert", "LINER")
                         else "radio_origin_unresolved" if r.get("radio_detection") in ("I", "unI", "I+unI")
                         else "not_detected_in_this_survey")
        out = {
            "galaxy": name, "activity_name": r.get("name", ""),
            "catalog_coverage": "matched" if match else "no_catalog_match",
            "match_method": "canonical_name" if key == match else "SIMBAD_alias" if match else "",
            "optical_class": activity_class, "radio_loud_flag": r.get("radio_loud_flag", ""),
            "jet_evidence": jet_class, "jet_evidence_source": jet.get("source_url", ""),
            "radio_detection": r.get("radio_detection", ""), "radio_morphology": r.get("radio_morphology", ""),
            "radio_core_limit": r.get("log_core_nuLnu_erg_s_limit", ""),
            "log_core_Lnu_1p5_W_Hz_published": r.get("log_core_Lnu_W_Hz", ""),
            "log_core_Lnu_1p5_W_Hz_sparc_distance": "", "radio_core_power_decade": "",
            "radio_core_rank_among_detections": "", "radio_core_rank_among_AGN_detections": "",
            "radio_distance_mpc": r.get("distance_mpc", ""),
            "log_oiii_erg_s": r.get("log_oiii_erg_s", ""), "oiii_limit": r.get("log_oiii_erg_s_limit", ""),
            "log_eddington_ratio": r.get("log_eddington_ratio", ""), "eddington_limit": r.get("log_eddington_ratio_limit", ""),
            "log_xray_03_10_erg_s_published": x.get("log_L_03_10_erg_s", ""),
            "log_xray_03_10_erg_s_sparc_distance": "", "xray_limit": x.get("log_L_03_10_erg_s_limit", ""),
            "xray_error_dex": x.get("log_L_03_10_erg_s_error_dex", ""),
            "xray_distance_mpc": x.get("distance_mpc", ""), "chandra_obsid": x.get("chandra_obsid", ""),
            "sparc_quality": meta["quality"], "inclination_deg": meta["inclination_deg"],
            "distance_mpc": meta["distance_mpc"], "disk_scale_kpc": meta["disk_scale_kpc"],
            "luminosity_1e9_lsun": meta["luminosity_1e9_lsun"], "hi_mass_1e9_msun": meta["hi_mass_1e9_msun"],
            "points": len(curve), "r_min_kpc": curve[0, 0], "r_max_kpc": curve[-1, 0],
            "rotation_quality_screen": int(meta["quality"] <= 2 and meta["inclination_deg"] >= 30),
            "fit_status": f.get("status", "not_in_selected_run"), "k": f.get("k", ""),
            "rms_kms": f.get("rms_kms", ""), "chi2": f.get("chi2", ""),
            "condition_number": f.get("condition_number", ""),
            "coefficient_bound_reached": f.get("coefficient_bound_reached", ""),
            **{k: "" for k in metric_names},
        }
        if match:
            value = float(r["log_core_Lnu_W_Hz"]) + 2*math.log10(meta["distance_mpc"]/float(r["distance_mpc"]))
            out["log_core_Lnu_1p5_W_Hz_sparc_distance"] = value
            if core_detection:
                out["radio_core_power_decade"] = f"10^{math.floor(value)} to 10^{math.floor(value)+1} W/Hz"
            if out["log_xray_03_10_erg_s_published"] != "":
                out["log_xray_03_10_erg_s_sparc_distance"] = float(out["log_xray_03_10_erg_s_published"]) + 2*math.log10(meta["distance_mpc"]/float(x["distance_mpc"]))
            baryons = baryonic_components(curve, ml_disk, ml_bulge)
            for c, b in zip(curve, baryons):
                points.append(dict(galaxy=name, R_kpc=c[0], Vobs_kms=c[1], eV_kms=c[2],
                                   v2_gas_kms2=b[0], v2_disk_kms2=b[1], v2_bulge_kms2=b[2],
                                   v2_baryons_kms2=b.sum()))
        if f.get("status") == "complete":
            if not math.isclose(float(f["distance_mpc"]), meta["distance_mpc"]):
                raise ValueError(f"{name}: fit and SPARC distances disagree")
            out.update(concentration(weights[name]))
            absolute_mass = math.fsum(abs(float(w["signed_mass_msun"])) for w in weights[name])
            for w in weights[name]:
                joined_weights.append({"galaxy": name, **{k:v for k,v in w.items() if k != "galaxy"},
                                       "abs_mass_fraction": abs(float(w["signed_mass_msun"]))/absolute_mass if absolute_mass else "",
                                       "scale_over_disk_scale": float(w["a0_star_kpc"])/meta["disk_scale_kpc"],
                                       "activity_name": out["activity_name"], "optical_class": activity_class,
                                       "jet_evidence": jet_class, "radio_core_limit": out["radio_core_limit"]})
        rows.append(out)
    # Only detections receive ranks. Limits and missing observations remain distinct.
    detections = [r for r in rows if r["radio_core_limit"] == "detection"]
    for label, subset in (("radio_core_rank_among_detections", detections),
                          ("radio_core_rank_among_AGN_detections", [r for r in detections if r["optical_class"] in ("Seyfert", "LINER")])):
        for rank, row in enumerate(sorted(subset, key=lambda r: -r["log_core_Lnu_1p5_W_Hz_sparc_distance"]), 1):
            row[label] = rank
    return rows, joined_weights, points, {"sources": manifest, "fit_run": run,
        "fit_input_sha256": fit_inputs, "ml_disk": ml_disk, "ml_bulge": ml_bulge}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "reports/agn")
    parser.add_argument("--fits-dir", type=Path, help="SPARC run directory with run.json, summary.csv and per-galaxy fit_weights.csv; default: bundled reference fits")
    args = parser.parse_args(argv)
    out = args.out_dir.resolve()
    if out == ROOT or DATA.resolve() == out or DATA.resolve() in out.parents:
        parser.error("output must be separate from project/source data")
    rows, weights, points, info = build_catalog(args.fits_dir)
    matched = [r for r in rows if r["catalog_coverage"] == "matched"]
    counts = {
        "sparc_total": len(rows), "catalog_matches": len(matched),
        "alias_matches": sum(r["match_method"] == "SIMBAD_alias" for r in rows),
        "optical_AGN": sum(r["optical_class"] in ("Seyfert", "LINER") for r in matched),
        "radio_core_detections": sum(r["radio_core_limit"] == "detection" for r in matched),
        "radio_core_upper_limits": sum(r["radio_core_limit"] == "upper" for r in matched),
        "AGN_radio_core_detections": sum(r["radio_core_rank_among_AGN_detections"] != "" for r in matched),
        "resolved_jet_lobe_priorities": sum(r["jet_evidence"] == "resolved_jet_lobes" for r in rows),
        "xray_detections": sum(r["xray_limit"] == "detection" for r in matched),
        "xray_upper_limits": sum(r["xray_limit"] == "upper" for r in matched),
        "xray_unobserved": sum(r["xray_limit"] == "" for r in matched),
        "rotation_points": len(points),
    }
    info["counts"] = counts
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "sparc175_activity_weights.csv", rows)
    write_csv(out / "matched51_activity_weights.csv", matched)
    write_csv(out / "orbital_weights.csv", weights)
    write_csv(out / "rotation_points.csv", points)
    atomic_write(out / "galaxies.txt", lambda f: f.write(",".join(r["galaxy"] for r in matched)+"\n"))
    atomic_write(out / "provenance.json", lambda f: json.dump(info, f, indent=2, allow_nan=False))
    print(json.dumps(counts, indent=2))
    print(f"Wrote {out}; strength ranks use measured radio cores only. Upper limits are not zero.")


if __name__ == "__main__":
    main()
