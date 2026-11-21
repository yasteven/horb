#!/usr/bin/env python3
"""Render an audited six-model comparison and its matching signed-density map."""
import argparse
import copy
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
import numpy as np

from fit_signed_density import atomic_write, fingerprint, write_csv
from horb_signed import load_basis, load_observations, baryons_function
from horb_highlights import MODELS, LABELS, fit_baseline, metrics, model_v2, load_density

ROOT = Path(__file__).resolve().parents[1]
GALAXIES = {
    "mw": {"title": "Milky Way", "source": "Sofue (2020)", "number": "01",
           "rc": "data/milky_way/sofue_2020_standard_rc.csv",
           "baryons": "data/milky_way/sofue_literature_baryons.csv",
           "note": "Fixed Sofue bulge + sphericalized exponential disk",
           "sources": ["https://doi.org/10.3390/galaxies8020037", "docs/signed_density.md"]},
    "m31": {"title": "Andromeda", "source": "Corbelli et al. (2010)", "number": "02",
            "rc": "reports/m31_inputs/rotation_curve.csv", "baryons": "reports/m31_inputs/baryons.csv",
            "note": "Fixed Geehan bulge + thin disk + central mass; no separate gas",
            "sources": ["https://arxiv.org/abs/0912.4133", "https://arxiv.org/abs/astro-ph/0501240"]},
}
BG, PANEL, INK, MUTED, GRID = "#0b1220", "#101b2d", "#f0f5fc", "#9cacc2", "#263448"
COLORS = dict(zip(MODELS, ["#5de4d0", "#c3a2ff", "#ffb65e", "#75a8ff", "#f487b4", "#a4b0c1"]))
STYLES = dict(zip(MODELS, ["-", "-", (0, (5, 2)), (0, (6, 2, 1, 2)), (0, (1, 1.5)), (0, (4, 3))]))


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "figure.facecolor": BG,
                         "axes.facecolor": PANEL, "axes.edgecolor": GRID, "text.color": INK,
                         "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                         "grid.color": GRID, "grid.alpha": .65, "axes.titleweight": "bold",
                         "savefig.facecolor": BG, "svg.fonttype": "none", "pdf.fonttype": 42})


def header(fig, galaxy, category, subtitle):
    g = GALAXIES[galaxy]
    fig.text(.07, .94, f"HORB   /   {category.upper()}", color=COLORS["superposition"], size=11, weight="bold")
    fig.text(.95, .94, g["number"], ha="right", color=MUTED, size=11)
    fig.text(.07, .872, g["title"], size=31, weight="bold")
    fig.text(.07, .828, subtitle, color=MUTED, size=11)


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    # Vector copies are useful for papers, slides and zooming; keep them out of git.
    pdf = ROOT / "plots" / (path.stem + ".pdf")
    pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf, metadata={"Creator": "HORB", "CreationDate": None, "ModDate": None})
    plt.close(fig)
    print(f"wrote {path} and {pdf}")


def rotation_figure(galaxy, r, curves, obs, predictions, scores, report, output):
    style()
    fig = plt.figure(figsize=(14, 8.2))
    g = GALAXIES[galaxy]
    lo, hi = report["r_fit_requested_kpc"]
    header(fig, galaxy, "Rotation curve", f"{g['source']}  ·  {len(obs)} observations  ·  {lo:g}–{hi:g} kpc  ·  shared fixed baryons")
    ax = fig.add_axes([.07, .355, .64, .425])
    res = fig.add_axes([.07, .17, .64, .145], sharex=ax)
    # Draw the selected superposition last so its path stays visible at overlaps.
    for name in reversed(MODELS):
        ax.plot(r, curves[name], color=COLORS[name], ls=STYLES[name],
                lw=3 if name == "superposition" else 1.8, zorder=3 if name == "superposition" else 2)
        if name != "baryons":
            res.plot(obs["R_kpc"], predictions[name]-obs["Vobs_kms"], color=COLORS[name],
                     ls=STYLES[name], lw=2.3 if name == "superposition" else 1.2)
    ax.errorbar(obs["R_kpc"], obs["Vobs_kms"], yerr=obs["eV_kms"], fmt="o", ms=4,
                mfc=INK, mec=BG, mew=.7, ecolor="#b9c4d3", elinewidth=.8, capsize=2, zorder=5)
    res.axhline(0, color=INK, alpha=.5, lw=.8)
    res.fill_between(obs["R_kpc"], -obs["eV_kms"], obs["eV_kms"], color=INK, alpha=.055)
    for a in (ax, res):
        a.grid(axis="y")
        a.spines[["top", "right"]].set_visible(False)
        a.set_xlim(lo, hi)
    ax.tick_params(labelbottom=False)
    ax.set_ylabel("Circular speed  [km/s]", labelpad=12)
    ax.set_ylim(0, max(np.max(obs["Vobs_kms"]+obs["eV_kms"]), max(np.max(v) for v in curves.values()))*1.08)
    res.set(xlabel="Galactocentric radius  [kpc]", ylabel="Δv  [km/s]")
    res.yaxis.labelpad = 12
    # Residuals show all five halo models; baryons-only is in the main panel.
    limit = max(10, max(np.max(np.abs(predictions[n]-obs["Vobs_kms"])) for n in MODELS if n != "baryons")*1.15)
    res.set_ylim(-limit, limit)
    fig.text(.755, .777, "MODEL", color=MUTED, size=9, weight="bold")
    fig.text(.951, .777, "RMS · km/s", color=MUTED, size=9, ha="right")
    for i, name in enumerate(MODELS):
        y = .727 - i*.077
        fig.add_artist(plt.Line2D([.755, .78], [y, y], transform=fig.transFigure,
                                 color=COLORS[name], ls=STYLES[name], lw=2.5))
        fig.text(.79, y, LABELS[name], color=COLORS[name], size=10, va="center")
        fig.text(.951, y-.025, f"{scores[name]['rms_kms']:.2f}", ha="right", size=12, weight="bold")
        description = {"superposition": f"{report['k']} signed-weight terms", "single": r"$3d_{z^2}$ · fitted mass + scale",
                       "nfw": "fitted mass + scale", "burkert": "fitted mass + scale",
                       "einasto": "fitted mass + scale + shape", "baryons": "fixed literature model"}[name]
        if name == "einasto" and scores[name]["bound_reached"]:
            description = f"α = {scores[name]['parameters'][2]:g} (parameter bound)"
        fig.text(.79, y-.026, description, color=MUTED, size=8)
    fig.text(.755, .19, "●  Observations with published errors\n\nΔv = model − observations\nShading = published error bars", size=9, color=MUTED, linespacing=1.4)
    fig.text(.07, .080, g["note"], color=MUTED, size=9)
    fig.text(.07, .049, "In-sample fits. Superposition selected by weighted χ² over the documented ≤4-term dictionary; model flexibility differs.", color=MUTED, size=8.7)
    save(fig, output)


def density_figure(galaxy, x, z, rho, report, output):
    style()
    fig = plt.figure(figsize=(14, 8.2))
    header(fig, galaxy, "Fitted dark density", "The exact superposition from the rotation plot above  ·  x–z slice through y = 0")
    # Signed log color: negative cyan, zero dark, positive amber. No abs(rho).
    cmap = LinearSegmentedColormap.from_list("horb_signed", ["#66e4f0", "#146d94", PANEL, "#b94f31", "#ffde83"])
    limit = 10**np.ceil(np.log10(max(float(np.max(np.abs(rho))), 1e-20)))
    threshold = limit*1e-5
    norm = SymLogNorm(linthresh=threshold, linscale=.7, vmin=-limit, vmax=limit, base=10)
    ax = fig.add_axes([.07, .19, .50, .59])
    im = ax.pcolormesh(x, z, rho, norm=norm, cmap=cmap, shading="auto", rasterized=True)
    ax.set(xlabel="x  [kpc]", ylabel="z  [kpc]", aspect="equal")
    ax.axhline(0, color=INK, lw=.6, alpha=.25)
    ax.axvline(0, color=INK, lw=.6, alpha=.25)
    zoom = min(12., max(abs(x[0]), abs(x[-1]))/3)
    ax.add_patch(plt.Rectangle((-zoom, -zoom), 2*zoom, 2*zoom, fill=False, edgecolor=INK, alpha=.7, lw=.8, ls="--"))
    cbax = fig.add_axes([.123, .11, .394, .022])
    cb = fig.colorbar(im, cax=cbax, orientation="horizontal", ticks=[-limit, -limit/100, 0, limit/100, limit])
    cb.ax.tick_params(labelsize=9)
    cb.set_label(r"Signed density  [$M_\odot\,/\,\mathrm{kpc}^3$]  ·  symmetric log", size=9)
    inset = fig.add_axes([.665, .455, .245, .325])
    inset.pcolormesh(x, z, rho, norm=norm, cmap=cmap, shading="auto", rasterized=True)
    inset.set(xlim=(-zoom, zoom), ylim=(-zoom, zoom), aspect="equal", title="Central detail · same color scale")
    inset.set_xlabel("x [kpc]", fontsize=9)
    inset.set_ylabel("z [kpc]", fontsize=9)
    inset.tick_params(labelsize=9)
    inset.title.set_fontsize(10)
    fig.text(.635, .365, r"$\rho_{\rm DM}=\sum_i w_i\,M_{\rm ref}\,|R_{n\ell}Y_{\ell m}|^2$", size=17)
    fig.text(.635, .315, "ORBITAL", color=MUTED, size=9, weight="bold")
    fig.text(.79, .315, "a₀ · kpc", color=MUTED, size=9, ha="right")
    fig.text(.945, .315, "SIGNED WEIGHT", color=MUTED, size=9, ha="right")
    for i, term in enumerate(report["components"]):
        y = .277-i*.035
        fig.text(.635, y, term["state"].replace("_", " "), size=11)
        fig.text(.79, y, f"{term['a0_star_kpc']:g}", size=11, ha="right")
        fig.text(.945, y, f"{term['weight']:+.6f}", size=11, ha="right", color="#ffde83" if term["weight"] >= 0 else "#66e4f0")
    fig.text(.635, .100, r"$M_{\rm ref}=10^{11}\,M_\odot$" + "  ·  amber + / cyan −\nWeights are not normalized probabilities.", size=10, color=MUTED, linespacing=1.7)
    fig.text(.07, .027, "Fitted model density; rotation data constrain disk-plane force and do not uniquely determine a 3D mass distribution.", color=MUTED, size=9)
    save(fig, output)


def validate_weights(path, terms):
    with Path(path).open() as f:
        rows = list(csv.DictReader(f))
    if len(rows) != len(terms):
        raise ValueError("weights CSV and fit report have different component counts")
    for row, term in zip(rows, terms):
        for key in ("basis_id", "n", "l", "m", "a0_star_kpc", "ref_mass_msun", "weight"):
            if not np.isclose(float(row[key]), term[key], rtol=1e-13, atol=1e-14):
                raise ValueError(f"weights CSV does not match fit report: {key}")


def update_readme_scores():
    """Keep the displayed table tied to the newly rendered figures."""
    readme = ROOT/"README.md"
    text = readme.read_text()
    start, end = "<!-- HORB_HIGHLIGHT_SCORES_START -->", "<!-- HORB_HIGHLIGHT_SCORES_END -->"
    if text.count(start) != 1 or text.count(end) != 1:
        raise ValueError("README score markers missing or duplicated; cannot update its score table")
    records = [json.loads((ROOT/f"data/highlights/{n}_comparison.json").read_text()) for n in ("mw", "m31")]
    fits = [json.loads((ROOT/f"data/highlights/{n}_fit.json").read_text()) for n in ("mw", "m31")]
    columns = [f"{d['galaxy']} · {d['fit_range_kpc'][0]:g}–{d['fit_range_kpc'][1]:g} kpc" for d in records]
    lines = ["| In-sample RMS (km/s) | " + " | ".join(columns) + " |", "|---|---:|---:|"]
    labels = [f"HORB superposition · {fits[0]['k']}/{fits[1]['k']} terms (MW/M31)",
              "HORB single `3d_z2` · fitted mass + scale", "NFW · fitted mass + scale", "Burkert · fitted mass + scale",
              "Einasto · fitted mass + scale + shape", "Baryons alone · fixed"]
    if fits[0]["k"] == fits[1]["k"]:
        labels[0] = f"HORB superposition · {fits[0]['k']} terms"
    for name, label in zip(MODELS, labels):
        values = [f"{d['scores'][name]['rms_kms']:.2f}" for d in records]
        if name == "superposition": values = [f"**{v}**" for v in values]
        lines.append("| " + " | ".join([label, *values]) + " |")
    before, rest = text.split(start)
    _, after = rest.split(end)
    atomic_write(readme, lambda f: f.write(before+start+"\n\n"+"\n".join(lines)+"\n\n"+end+after))


def prepare(galaxy, report, basis_path):
    """Recompute every plotted speed and score, refusing stale/mismatched fits."""
    g = GALAXIES[galaxy]
    basis = load_basis(basis_path)
    _, obs = load_observations(ROOT/g["rc"], *report["r_fit_requested_kpc"], basis.radii)
    baryons = baryons_function(ROOT/g["baryons"])
    lookup = {m["basis_id"]: i for i, m in enumerate(basis.metadata)}
    indices, weights = [], []
    for term in report["components"]:
        index = lookup[term["basis_id"]]
        for key, value in basis.metadata[index].items():
            if term[key] != value:
                raise ValueError(f"selected component {term['basis_id']} has changed {key}")
        indices.append(index)
        weights.append(term["weight"])
    r_obs = obs["R_kpc"]
    b_obs = baryons(r_obs)
    dark = basis.at(r_obs)[:, indices] @ weights
    snapshot = np.column_stack([r_obs, obs["Vobs_kms"], obs["eV_kms"], b_obs, dark])
    if "fitted_observations" in report:
        if not np.allclose(snapshot, report["fitted_observations"], rtol=1e-9, atol=1e-8):
            raise ValueError("fit inputs or selected forces have changed; use --refit")
    else:
        for key, path in [("observations", g["rc"]), ("baryons", g["baryons"])]:
            if fingerprint(ROOT/path)["sha256"] != report["inputs"][key]["sha256"]:
                raise ValueError(f"fit input {key} has changed; refit before plotting")
    # Recheck the same dense positivity guard, including nodes between the samples.
    guard = np.unique(np.concatenate([basis.radii, r_obs, *[
        basis.radii[:-1]*(1-t)+basis.radii[1:]*t for t in (.25, .5, .75)]]))
    if np.any(baryons(guard)+basis.at(guard)[:, indices] @ weights < 0):
        raise ValueError("selected superposition has negative total circular force")
    r = np.unique(np.concatenate([np.linspace(*report["r_fit_requested_kpc"], 1000), r_obs]))
    predictions = {"baryons": np.sqrt(b_obs), "superposition": np.sqrt(b_obs+dark)}
    curves = {"baryons": np.sqrt(baryons(r)), "superposition": np.sqrt(baryons(r)+basis.at(r)[:, indices] @ weights)}
    scores = {key: metrics(value, obs["Vobs_kms"], obs["eV_kms"]) for key, value in predictions.items()}
    for key in ("rms_kms", "chi2"):
        if not np.isclose(scores["superposition"][key], report[key], rtol=1e-7, atol=1e-9):
            raise ValueError(f"recomputed {key} does not match selected fit")
    for name in ("single", "nfw", "burkert", "einasto"):
        result = fit_baseline(name, r_obs, obs["Vobs_kms"], obs["eV_kms"], b_obs, report["r_fit_requested_kpc"][1])
        scores[name] = result
        for radii, b, dest in [(r_obs, b_obs, predictions), (r, baryons(r), curves)]:
            dest[name] = np.sqrt(b+model_v2(name, radii, result["log10_parameters"], result["reference_radius_kpc"]))
    return r, curves, obs, predictions, scores, snapshot


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--galaxy", choices=GALAXIES, required=True)
    p.add_argument("--fit-report", type=Path, required=True)
    p.add_argument("--weights", type=Path, required=True)
    p.add_argument("--basis-csv", type=Path, required=True)
    p.add_argument("--density-csv", type=Path, required=True)
    args = p.parse_args()
    report = json.loads(args.fit_report.read_text())
    validate_weights(args.weights, report["components"])
    r, curves, obs, predictions, scores, snapshot = prepare(args.galaxy, report, args.basis_csv)
    x, z, rho = load_density(args.density_csv, report["components"])
    name = args.galaxy
    rotation_figure(name, r, curves, obs, predictions, scores, report, ROOT/f"images/{name}_rotation_highlight.png")
    density_figure(name, x, z, rho, report, ROOT/f"images/{name}_density_highlight.png")
    record = ROOT/"data/highlights"
    record.mkdir(parents=True, exist_ok=True)
    portable = copy.deepcopy(report)
    portable["fitted_observations_columns"] = ["R_kpc", "Vobs_kms", "eV_kms", "v2_baryons_kms2", "v2_dark_signed_kms2"]
    portable["fitted_observations"] = snapshot.tolist()
    for item in portable["inputs"].values():
        item["path"] = str(Path(item["path"]).relative_to(ROOT)) if Path(item["path"]).is_relative_to(ROOT) else item["path"]
    portable["config"]["basis_csv"] = f"reports/{name}_showcase_basis.csv"
    portable["inputs"]["basis"]["path"] = portable["config"]["basis_csv"]
    atomic_write(record/f"{name}_fit.json", lambda f: json.dump(portable, f, indent=2, allow_nan=False))
    weight_fields = ["basis_id", "state", "n", "l", "m", "a0_star_kpc", "ref_mass_msun", "weight", "signed_mass_msun"]
    write_csv(record/f"{name}_weights.csv", weight_fields, [[t[k] for k in weight_fields] for t in report["components"]])
    comparison = {"galaxy": GALAXIES[name]["title"], "sources": GALAXIES[name]["sources"], "scores": scores,
                  "fit_range_kpc": report["r_fit_requested_kpc"], "baryons": GALAXIES[name]["note"],
                  "selection": "Lowest converged weighted chi2 in the enumerated <=4-term dictionary; no model-probability claim.",
                  "density": {"slice": "y=0", "extent_kpc": float(x[-1]), "bins": len(x),
                              "min_msun_per_kpc3": float(rho.min()), "max_msun_per_kpc3": float(rho.max())},
                  "validation": "Weights, signed-density CSV, input snapshot, fitted scores and total-force guard checked before rendering."}
    atomic_write(record/f"{name}_comparison.json", lambda f: json.dump(comparison, f, indent=2, allow_nan=False))
    write_csv(record/f"{name}_scores.csv", ["model", "n", "rms_kms", "chi2", "chi2_per_point"],
              [[LABELS[m], scores[m]["n"], scores[m]["rms_kms"], scores[m]["chi2"], scores[m]["chi2"]/scores[m]["n"]] for m in MODELS])
    write_csv(record/f"{name}_curves.csv", ["r_kpc", *[f"v_{m}_kms" for m in MODELS]], np.column_stack([r, *[curves[m] for m in MODELS]]))
    update_readme_scores()
    print(json.dumps({m: round(scores[m]["rms_kms"], 4) for m in MODELS}))


if __name__ == "__main__":
    main()
