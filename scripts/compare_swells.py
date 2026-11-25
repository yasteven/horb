#!/usr/bin/env python3
"""Compare positive-only and signed HORB selections from saved rotation fits."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from fit_swells import (ROOT, DATA, MODELS, LABELS, COLORS, fingerprints,
                        force_grid, total_force, write_json)
from horb_highlights import metrics
from horb_swells.models import observing_matrix
from horb_swells.lensing import orbital_surface_modes, FourierLens, CombinedLens


def component_key(components):
    return json.dumps(components, sort_keys=True, allow_nan=False)


def load_family(directory, nonnegative):
    report = json.loads((directory/"rotation_fit.json").read_text())
    lensing = json.loads((directory/"lensing_predictions.json").read_text())
    if (report["fit_inputs_sha256"] != fingerprints()
            or lensing["fit_inputs_sha256"] != report["fit_inputs_sha256"]
            or report["lensing_used_in_fit"] or lensing["lensing_used_in_fit"]
            or report["search"]["nonnegative"] != nonnegative):
        raise ValueError(f"{directory}: incompatible inputs, sign constraint or fitting provenance")
    if lensing["source"] != json.loads((DATA/"j2141_lensing.json").read_text()):
        raise ValueError(f"{directory}: lensing summary changed; regenerate predictions")
    return report, lensing


def candidate_pool(report, lensing):
    """Use the recorded pool, including the rotation winner, without refitting."""
    best = report["models"]["superposition"]
    fits = {component_key(c["components"]): c
            for c in [best, *report["alternative_rotation_fits"]]}
    predictions = [{**lensing["models"]["superposition"], "components": best["components"]},
                   *lensing["alternative_rotation_fits"]]
    pool, seen = [], set()
    for prediction in predictions:
        key = component_key(prediction["components"])
        if key not in fits:
            raise ValueError("lensing candidate has no matching saved rotation fit")
        if key in seen:
            continue
        seen.add(key)
        fit = fits[key]
        if report["search"]["nonnegative"] and any(t["weight"] < -1e-10 for t in fit["components"]):
            raise ValueError("negative coefficient in positive-only candidate pool")
        if not np.isfinite([fit["chi2"], fit["rms_kms"]]).all():
            raise ValueError("nonfinite rotation scores")
        theta = prediction["effective_einstein_radius_arcsec"]
        if theta is not None and (not np.isfinite(theta) or theta <= 0):
            raise ValueError("invalid Einstein radius")
        pool.append({"rotation_chi2": fit["chi2"], "rotation_rms_kms": fit["rms_kms"],
                     "n_rotation": fit["n"], "components": fit["components"],
                     "effective_einstein_radius_arcsec": theta,
                     "maximum_outward_crossings_per_ray": prediction["maximum_outward_crossings_per_ray"]})
    return pool


def select_family(pool, observed, family):
    rotation = min(pool, key=lambda c: c["rotation_chi2"])
    valid = [c for c in pool if c["effective_einstein_radius_arcsec"] is not None
             and c["maximum_outward_crossings_per_ray"] == 1]
    if not valid:
        raise ValueError(f"{family}: no single-branch lensing radius among saved candidates")
    lensing = min(valid, key=lambda c: (abs(c["effective_einstein_radius_arcsec"]-observed),
                                       c["rotation_chi2"]))
    return [{**candidate, "family": family, "selection": selection,
             "selection_uses_lensing": selection == "closest_lensing"}
            for selection, candidate in (("best_rotation", rotation), ("closest_lensing", lensing))]


def compare(signed, positive):
    sr, sl = signed
    pr, pl = positive
    for key in ("fit_inputs_sha256", "photometry", "distances", "baryons", "observations", "fit_n",
                "stellar_mass_offset_dex"):
        if sr[key] != pr[key]:
            raise ValueError(f"positive and signed runs have different {key}; use matched assumptions")
    if sl["source"] != pl["source"]:
        raise ValueError("positive and signed runs use different lensing summaries")
    observed = sl["source"]["circularized_einstein_radius_arcsec"]
    rows, pools = [], {}
    for family, report, lensing in (("positive_only", pr, pl), ("signed", sr, sl)):
        pool = candidate_pool(report, lensing)
        pools[family] = {"saved_candidates": len(pool),
                         "subsets_evaluated_in_rotation_search": report["search"]["subsets_evaluated"],
                         "candidate_scope": lensing["alternatives_note"]}
        rows.extend(select_family(pool, observed, family))
    for name in MODELS[1:]:
        fit = sr["models"][name]
        rows.append({"family": name, "selection": "best_rotation", "selection_uses_lensing": False,
                     "rotation_rms_kms": fit["rms_kms"], "rotation_chi2": fit["chi2"],
                     "n_rotation": fit["n"],
                     "effective_einstein_radius_arcsec": sl["models"][name]["effective_einstein_radius_arcsec"]})
    return {"schema_version": 1, "galaxy": sr["galaxy"], "source": sl["source"],
            "fit_inputs_sha256": sr["fit_inputs_sha256"], "pools": pools,
            "rotation_selection_statistic": "minimum weighted rotation chi2",
            "lensing_selection_statistic": "minimum absolute Einstein-radius error among saved rotation-fitted candidates",
            "selection_scope": "Best within each recorded pool, not a global lensing optimization. Weights remain rotation-fitted; lens-selected rows use the observed radius for candidate selection.",
            "density_map_selection": {"family": "signed", "selection": "closest_lensing"},
            "rows": rows}


def audit_and_map(result, signed, positive, out):
    """Independently reconstruct every displayed HORB score and selected map."""
    radial = np.geomspace(1e-5, 10000, 3200)
    cache, selected_map = {}, None
    for row in result["rows"][:4]:
        report = positive[0] if row["family"] == "positive_only" else signed[0]
        meta, distance = report["photometry"], report["distances"]
        sigma = distance["sigma_critical_msun_kpc2"]
        observations = report["observations"]
        mask = np.array([o["use_for_fit"] for o in observations], dtype=bool)
        grid = force_grid()
        P = observing_matrix(np.array([o["r_arcsec"] for o in observations]), grid, meta)
        candidate = {**report, "models": {"superposition": {"components": row["components"]}}}
        prediction = P@np.sqrt(total_force("superposition", grid, candidate))
        score = metrics(prediction[mask], np.array([o["v_los_kms"] for o in observations])[mask],
                        np.array([o["error_kms"] for o in observations])[mask])
        if abs(score["chi2"]-row["rotation_chi2"]) > 1e-6:
            raise ValueError("selected candidate does not reproduce its saved rotation score")
        row["predicted_los_kms"] = prediction.tolist()
        surface = np.zeros((3, len(radial)))
        for term in row["components"]:
            key = tuple(term[k] for k in ("n", "l", "a0_star_kpc", "ref_mass_msun"))
            if key not in cache:
                cache[key] = orbital_surface_modes(radial, term, meta["cos_inclination"])
            surface += term["weight"]*cache[key]
        lens = CombinedLens(FourierLens(radial, surface/sigma), report["baryons"], sigma)
        critical = lens.critical_curve(distance["kpc_per_arcsec"], angles=256)
        value, recorded = critical["effective_einstein_radius_arcsec"], row["effective_einstein_radius_arcsec"]
        if ((value is None) != (recorded is None)
                or value is not None and abs(value-recorded) > 2e-4):
            raise ValueError("selected candidate lensing radius failed reconstruction")
        row["refined_einstein_radius_arcsec"] = value
        if row["family"] == "signed" and row["selection"] == "closest_lensing":
            axis = np.linspace(-3.5, 3.5, 241)
            xx, yy = np.meshgrid(axis*distance["kpc_per_arcsec"], axis*distance["kpc_per_arcsec"])
            dark = lens.dark.field(xx, yy)
            selected_map = {"x_arcsec": axis, "y_arcsec": axis, "kappa_dark": (dark[2]+dark[4])/2,
                            "phi_radians": np.array(critical["phi_radians"]),
                            "critical_radius_arcsec": np.array(critical["radius_kpc"])/distance["kpc_per_arcsec"]}
    np.savez_compressed(out/"j2141_selected_lens_map.npz", **selected_map)
    return selected_map


def plot_comparison(result, lens_map, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib import patheffects
    bg, panel, fg, muted = "#0b1320", "#101e30", "#edf3fc", "#9dadc3"
    family_colors = {**COLORS, "positive_only": "#b8e878", "signed": COLORS["superposition"]}
    plt.rcParams.update({"figure.facecolor": bg, "axes.facecolor": panel, "axes.labelcolor": muted,
                         "text.color": fg, "xtick.color": muted, "ytick.color": muted,
                         "axes.edgecolor": "#35465a", "font.size": 11})
    fig = plt.figure(figsize=(17, 9.5))
    left = fig.add_axes([.065, .265, .30, .537])
    right = fig.add_axes([.64, .265, .275, .537])
    fig.text(.065, .945, "HORB  /  ROTATION + LENSING", color=COLORS["superposition"], size=14, weight="bold")
    fig.text(.065, .89, "SDSS J2141−0001", size=27, weight="bold")
    fig.text(.065, .85, "Positive-only and signed superpositions · separate rotation and lensing selections", color=muted, size=12)
    selected = result["rows"][3]
    dark = lens_map["kappa_dark"]
    limit = max(float(np.max(np.abs(dark))), .01)
    im = left.imshow(dark, origin="lower", extent=(-3.5, 3.5, -3.5, 3.5), cmap="RdBu_r",
                     norm=SymLogNorm(linthresh=max(limit*.01, 1e-4), vmin=-limit, vmax=limit))
    phi = np.r_[lens_map["phi_radians"], 2*np.pi]
    radius = np.r_[lens_map["critical_radius_arcsec"], lens_map["critical_radius_arcsec"][0]]
    line, = left.plot(radius*np.cos(phi), radius*np.sin(phi), color=fg, lw=2)
    line.set_path_effects([patheffects.Stroke(linewidth=4, foreground=bg), patheffects.Normal()])
    left.set(xlabel="Major-axis offset [arcsec]", ylabel="Minor-axis offset [arcsec]",
             title="Signed HORB · closest lensing")
    colorbar = fig.colorbar(im, cax=fig.add_axes([.376, .265, .012, .537]))
    colorbar.set_label(r"$\kappa_{\rm HORB}=\Sigma_{\rm HORB}/\Sigma_{\rm crit}$", labelpad=7)
    fig.text(.065, .18, f"Map: θE = {selected['effective_einstein_radius_arcsec']:.3f}″  ·  rotation RMS = {selected['rotation_rms_kms']:.3f} km/s", size=11)
    fig.text(.065, .153, "Projected dark density + total critical curve · GR thin lens", color=muted, size=10)
    source = result["source"]
    observed = source["circularized_einstein_radius_arcsec"]
    right.axvspan(observed-source["radius_error_minus_arcsec"], observed+source["radius_error_plus_arcsec"], color=fg, alpha=.12)
    right.axvline(observed, color=fg, ls="--", lw=1.2)
    labels, values = [], []
    rows = [row for row in result["rows"] if row["family"] != "single"]
    for j, row in enumerate(rows):
        family = row["family"]
        color = family_colors[family]
        lens_selected = row["selection_uses_lensing"]
        if family in ("positive_only", "signed"):
            labels.append(f"HORB · {'positive only' if family == 'positive_only' else 'signed (+/−)'}\n{'closest lensing' if lens_selected else 'best rotation'}")
        else:
            labels.append(LABELS[family])
        value = row["effective_einstein_radius_arcsec"]
        if value is not None:
            right.scatter(value, j, color=color, s=62, marker="D" if lens_selected else "o", zorder=5)
            right.annotate(f"{value:.3f}″", (value, j), xytext=(0, 9), textcoords="offset points", ha="center", color=color, size=10)
            values.append(value)
        else:
            right.text(.5, j, "multiple / missing branches", ha="center", color=color, size=8)
        right.text(1.10, j, f"{row['rotation_rms_kms']:.3f}", transform=right.get_yaxis_transform(),
                   ha="center", va="center", color=color, size=11, weight="bold")
    right.set_yticks(range(len(rows)), labels)
    right.tick_params(axis="y", labelsize=10.5, length=0, pad=12)
    for label, row in zip(right.get_yticklabels(), rows):
        label.set_color(family_colors[row["family"]])
    right.set(xlabel="Area-equivalent Einstein radius [arcsec]", ylim=(len(rows)-.4, -.7),
              xlim=(max(0., min(values)-.12), max(1.02, max(values)+.08)))
    right.grid(axis="x", alpha=.12)
    right.axhline(3.5, color=muted, alpha=.3, lw=.8)
    right.text(1.10, 1.025, "Rotation RMS\n[km/s]", transform=right.transAxes, ha="center", color=muted, size=10)
    fig.text(.64, .195, "●  Best rotation (weighted χ²)     ◆  Closest lensing radius", color=muted, size=10)
    counts = result["pools"]
    fig.text(.065, .105, f"Selections from saved rotation fits: {counts['positive_only']['saved_candidates']} positive-only, {counts['signed']['saved_candidates']} signed. Diamonds use the lensing measurement for selection.", color=muted, size=10.5)
    fig.text(.065, .072, f"White band: published SIE radius {observed:.2f} (+{source['radius_error_plus_arcsec']:.2f} / −{source['radius_error_minus_arcsec']:.2f}) arcsec. Same six rotation measurements and fixed SPS baryons.", color=muted, size=10.5)
    fig.text(.065, .04, "Closest lensing refers to this saved candidate pool; no continuous lensing refit or observed-arc pixel likelihood.", color=muted, size=10.5)
    fig.savefig(out/"j2141_lensing.png", dpi=180)
    fig.savefig(out/"j2141_lensing.pdf")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--signed-dir", type=Path, default=DATA/"reference/j2141")
    p.add_argument("--positive-dir", type=Path, default=DATA/"reference/j2141_nonnegative")
    p.add_argument("--out-dir", type=Path, default=ROOT/"reports/swells/j2141_comparison")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args(argv)
    signed, positive = load_family(args.signed_dir, False), load_family(args.positive_dir, True)
    result = compare(signed, positive)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    lens_map = audit_and_map(result, signed, positive, args.out_dir)
    write_json(args.out_dir/"selection.json", result)
    fields = ["family", "selection", "selection_uses_lensing", "rotation_rms_kms", "rotation_chi2",
              "n_rotation", "effective_einstein_radius_arcsec"]
    with (args.out_dir/"summary.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result["rows"])
    if not args.no_plots:
        plot_comparison(result, lens_map, args.out_dir)
    for row in result["rows"][:4]:
        print(f"{row['family']} / {row['selection']}: RMS {row['rotation_rms_kms']:.3f} km/s; θE {row['effective_einstein_radius_arcsec']}")
    print(f"Wrote four HORB selections and baseline comparisons to {args.out_dir}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        sys.exit(f"SWELLS comparison: {exc}")
