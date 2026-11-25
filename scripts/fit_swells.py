#!/usr/bin/env python3
"""J2141 rotation-only fits followed by frozen-model lensing predictions."""
import os
# Small dense least-squares problems are faster and reproducible with one thread.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from horb_highlights import metrics, signed_density_xz
from horb_swells.models import (distances, baryon_components, baryon_v2, dictionary,
                               basis_matrix, baseline_v2, observing_matrix)
from horb_swells.fitting import fit_classical, ProjectedSearch
from horb_swells.lensing import (orbital_surface_modes, spherical_surface,
                                FourierLens, CombinedLens)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/swells"
MODELS = ("superposition", "single", "nfw", "burkert", "einasto", "baryons")
LABELS = {"superposition": "HORB superposition", "single": "HORB single orbital",
          "nfw": "NFW", "burkert": "Burkert", "einasto": "Einasto", "baryons": "Baryons alone"}
COLORS = {"superposition": "#50e0cf", "single": "#bf9bff", "nfw": "#ffb451",
          "burkert": "#70a8ff", "einasto": "#ee82af", "baryons": "#9cabbc"}


def read_inputs():
    meta = json.loads((DATA/"j2141_photometry.json").read_text())
    obs = np.atleast_1d(np.genfromtxt(DATA/"j2141_rotation.csv", delimiter=",", names=True))
    if (len(obs) != 10 or not np.all(np.isfinite(obs.tolist()))
            or np.any(obs["error_kms"] <= 0) or np.any(np.diff(obs["r_arcsec"]) <= 0)
            or not np.array_equal(obs["use_for_fit"].astype(bool), obs["r_arcsec"] > meta["fit_min_arcsec"])):
        raise ValueError("invalid J2141 observations or fitting mask")
    return meta, obs


def fingerprints():
    return {name: hashlib.sha256((DATA/name).read_bytes()).hexdigest()
            for name in ("j2141_photometry.json", "j2141_rotation.csv")}


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+"\n")


def force_grid(count=800):
    return np.r_[0., np.geomspace(1e-4, 100, count), np.geomspace(101, 1e4, 80)]


def total_force(model, r, report):
    baryons = baryon_v2(r, report["baryons"])
    result = report["models"][model]
    if model == "superposition":
        if result["components"]:
            baryons += basis_matrix(r, result["components"])@np.array([t["weight"] for t in result["components"]])
    elif model != "baryons":
        baryons += baseline_v2(model, r, result["log10_parameters"])
    if np.min(baryons) < -1e-6:
        raise ValueError(f"{model}: negative total circular force on evaluation grid")
    return np.maximum(baryons, 0)


def fit_rotation(meta, obs, args):
    mask = obs["use_for_fit"].astype(bool)
    grid, terms = force_grid(), dictionary()
    baryons = baryon_components(meta, args.stellar_mass_offset_dex)
    b = baryon_v2(grid, baryons)
    P = observing_matrix(obs["r_arcsec"], grid, meta)
    A = basis_matrix(grid, terms)
    v, e = obs["v_los_kms"][mask], obs["error_kms"][mask]
    search = ProjectedSearch(A, b, P[mask], v, e, max_weight=args.max_abs_weight,
                             nonnegative=args.nonnegative,
                             validation_matrix=basis_matrix(np.unique(np.r_[force_grid(1600), force_grid(3200)]), terms),
                             validation_baryons=baryon_v2(np.unique(np.r_[force_grid(1600), force_grid(3200)]), baryons))
    print("Fitting J2141 rotation measurements (6 outer points; lensing held out)...", flush=True)
    candidates = search.run(args.max_components, args.beam_width, args.budget, args.exhaustive)
    best = candidates[0]
    def components(candidate):
        return [{**terms[i], "weight": float(w)} for i, w in zip(candidate["subset"], candidate["weights"])]
    models = {"superposition": {**best, "components": components(best),
                                "parameter_count": len(best["subset"]),
                                "discrete_selection": "orbital identities and scales searched"}}
    print(f"HORB: RMS {best['rms_kms']:.4f} km/s; searched {len(search.cache)} subsets", flush=True)
    for name in ("single", "nfw", "burkert", "einasto"):
        models[name] = fit_classical(name, grid, b, P[mask], v, e)
        print(f"{LABELS[name]}: RMS {models[name]['rms_kms']:.4f} km/s", flush=True)
    models["baryons"] = {**metrics(P[mask]@np.sqrt(b), v, e), "parameter_count": 0}
    return {"schema_version": 1, "galaxy": meta["galaxy"], "fit_inputs_sha256": fingerprints(),
            "lensing_used_in_fit": False, "selection_statistic": "velocity-space weighted chi2",
            "observed_velocity_convention": "LOS; compare after inclination, PSF, slit and bin integration",
            "assumptions": ["fixed Chabrier SPS stellar masses plus requested common mass offset",
                            "fixed photometric geometry; aligned axisymmetric orbital density",
                            "Gaussian 0.60 arcsec seeing; exponential emitting thin disk",
                            "cold gas omitted: no measured gas mass map available",
                            "signed density allowed; inward total equatorial force guarded to 10000 kpc"],
            "stellar_mass_offset_dex": args.stellar_mass_offset_dex,
            "photometry": meta, "distances": distances(meta), "baryons": baryons,
            "observations": [{k: float(row[k]) for k in obs.dtype.names} for row in obs],
            "fit_n": int(mask.sum()), "excluded_n": int((~mask).sum()),
            "search": {"method": "exhaustive" if args.exhaustive else "beam", "beam_width": args.beam_width,
                       "budget": args.budget, "max_components": args.max_components,
                       "max_abs_weight": args.max_abs_weight, "nonnegative": args.nonnegative,
                       "dictionary": terms, "subsets_evaluated": len(search.cache), "failed_subsets": search.failed},
            "models": models,
            "alternative_rotation_fits": [{**c, "components": components(c)} for c in candidates[1:21]],
            "uncertainty_note": "Point fits conditional on baryons and geometry. Alternative subsets are not posterior samples or confidence intervals."}


def rotation_outputs(report, meta, obs, out):
    # Double both force and integration resolution for a numerical audit, without
    # changing parameters or reselecting the model against lensing information.
    grid = force_grid(1600)
    P = observing_matrix(obs["r_arcsec"], grid, meta, step_arcsec=.02)
    mask = obs["use_for_fit"].astype(bool)
    table = np.column_stack([obs[k] for k in obs.dtype.names]+[
        P@np.sqrt(total_force(name, grid, report)) for name in MODELS])
    header = ",".join(list(obs.dtype.names)+[f"v_{name}_los_kms" for name in MODELS])
    np.savetxt(out/"rotation_predictions.csv", table, delimiter=",", header=header, comments="", fmt="%.12g")
    report["numerical_audit"] = {}
    coarse_grid = force_grid()
    coarse_P = observing_matrix(obs["r_arcsec"], coarse_grid, meta)
    for i, name in enumerate(MODELS):
        fine = table[:, len(obs.dtype.names)+i]
        coarse = coarse_P@np.sqrt(total_force(name, coarse_grid, report))
        difference = float(np.max(np.abs(fine[mask]-coarse[mask])))
        report["numerical_audit"][name] = {"max_fitted_velocity_change_kms": difference,
                                         "refined_scores": metrics(fine[mask], obs["v_los_kms"][mask], obs["error_kms"][mask])}
        if difference > .15:
            raise ValueError(f"{name}: observing integration not converged ({difference:.3f} km/s)")
    return table


def lens_predictions(report, meta, out):
    # Only here, after every rotation parameter has been fixed, read the
    # observed lensing summary. It is never passed to a mass-model optimizer.
    observed = json.loads((DATA/"j2141_lensing.json").read_text())
    distance = report["distances"]
    sigma = distance["sigma_critical_msun_kpc2"]
    radial = np.geomspace(1e-5, 10000, 3200)
    cache = {}
    def modes(components):
        result = np.zeros((3, len(radial)))
        for term in components:
            key = (term["n"], term["l"], term["a0_star_kpc"], term["ref_mass_msun"])
            if key not in cache:
                cache[key] = orbital_surface_modes(radial, term, meta["cos_inclination"])
            result += term["weight"]*cache[key]
        return result
    result = {"schema_version": 1, "source": observed, "models": {},
              "fit_inputs_sha256": report["fit_inputs_sha256"], "lensing_used_in_fit": False,
              "gravity_assumption": "GR weak-field thin lens; fitted signed density sources both metric potentials equally",
              "comparison_scope": "Published SIE summary only. No observed-image pixel likelihood, PSF/dust/source reconstruction, or model evidence computed.",
              "external_shear": 0.0, "external_convergence": 0.0,
              "distances": distance, "conditional_on_fixed_baryons_and_geometry": True}
    fields, curves, all_modes = {}, {}, []
    for name in MODELS:
        fit = report["models"][name]
        if name == "superposition":
            surface = modes(fit["components"])
        elif name == "single":
            mass, scale = fit["parameters"]
            surface = modes([{"n": 3, "l": 2, "a0_star_kpc": scale, "ref_mass_msun": mass, "weight": 1.}])
        elif name == "baryons":
            surface = np.zeros((3, len(radial)))
        else:
            surface = np.zeros((3, len(radial)))
            surface[0] = spherical_surface(radial, name, fit["log10_parameters"])
        all_modes.extend(surface)
        lens = CombinedLens(FourierLens(radial, surface/sigma), report["baryons"], sigma)
        critical = lens.critical_curve(distance["kpc_per_arcsec"])
        curves[name] = critical
        fields[name] = lens
        result["models"][name] = {k: v for k, v in critical.items() if k not in ("phi_radians", "radius_kpc")}
        print(f"{LABELS[name]}: predicted effective Einstein radius {critical['effective_einstein_radius_arcsec']} arcsec", flush=True)
    # A small documented collection of near-best searched subsets exposes the
    # extrapolation sensitivity. It is not selected by lensing agreement.
    alternatives = []
    best_chi = report["models"]["superposition"]["chi2"]
    for fit in report["alternative_rotation_fits"]:
        if fit["chi2"] > best_chi+1:
            continue
        lens = CombinedLens(FourierLens(radial, modes(fit["components"])/sigma), report["baryons"], sigma)
        curve = lens.critical_curve(distance["kpc_per_arcsec"], angles=128)
        alternatives.append({"chi2": fit["chi2"], "rms_kms": fit["rms_kms"], "components": fit["components"],
                             "effective_einstein_radius_arcsec": curve["effective_einstein_radius_arcsec"],
                             "maximum_outward_crossings_per_ray": curve["maximum_outward_crossings_per_ray"]})
    result["alternative_rotation_fits"] = alternatives
    result["alternatives_note"] = "Up to 20 next-best searched subsets with delta chi2 <= 1; illustrates degeneracy, not a confidence interval."
    np.savetxt(out/"surface_density_modes.csv", np.column_stack([radial, *all_modes]), delimiter=",", comments="",
               header="r_kpc,"+",".join(f"sigma_{name}_m{m}_msun_kpc2" for name in MODELS for m in (0, 2, 4)))
    phi = curves["superposition"]["phi_radians"]
    np.savetxt(out/"critical_curves.csv", np.column_stack([phi]+[curves[name]["radius_kpc"] for name in MODELS]),
               delimiter=",", header="phi_radians,"+",".join(f"r_{name}_kpc" for name in MODELS), comments="")
    axis = np.linspace(-3.5, 3.5, 241)
    xx, yy = np.meshgrid(axis*distance["kpc_per_arcsec"], axis*distance["kpc_per_arcsec"])
    f = fields["superposition"].field(xx, yy)
    dark = fields["superposition"].dark.field(xx, yy)
    np.savez_compressed(out/"j2141_lens_map.npz", x_arcsec=axis, y_arcsec=axis,
                        alpha_x_arcsec=f[0]/distance["kpc_per_arcsec"], alpha_y_arcsec=f[1]/distance["kpc_per_arcsec"],
                        kappa=(f[2]+f[4])/2, kappa_dark=(dark[2]+dark[4])/2,
                        gamma1=(f[2]-f[4])/2, gamma2=f[3],
                        inverse_magnification=(1-f[2])*(1-f[4])-f[3]**2)
    write_json(out/"lensing_predictions.json", result)
    with (out/"summary.csv").open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(["model", "rms_los_kms", "chi2_rotation", "n_rotation", "effective_einstein_radius_arcsec"])
        for name in MODELS:
            fit = report["models"][name]
            writer.writerow([name, fit["rms_kms"], fit["chi2"], fit["n"],
                             result["models"][name]["effective_einstein_radius_arcsec"]])
    return result, curves


def plot_results(report, lensing, curves, meta, obs, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    bg, panel, fg, muted = "#0b1320", "#101e30", "#edf3fc", "#9dadc3"
    plt.rcParams.update({"figure.facecolor": bg, "axes.facecolor": panel, "axes.labelcolor": muted,
                         "text.color": fg, "xtick.color": muted, "ytick.color": muted,
                         "axes.edgecolor": "#35465a", "font.size": 11})
    scale = report["distances"]["kpc_per_arcsec"]
    mask = obs["use_for_fit"].astype(bool)
    fig, ax = plt.subplots(figsize=(13.5, 7.5))
    fig.subplots_adjust(left=.08, right=.72, top=.81, bottom=.16)
    fig.text(.08, .93, "HORB  /  SWELLS", color=COLORS["superposition"], size=13, weight="bold")
    fig.text(.08, .87, "SDSS J2141−0001", size=25, weight="bold")
    fig.text(.08, .82, "Rotation fit · six outer measurements · lensing held out", color=muted)
    xx = np.linspace(0, obs["r_arcsec"].max(), 180)
    grid = force_grid(1600)
    P = observing_matrix(xx, grid, meta)
    ax.axvspan(0, meta["fit_min_arcsec"]*scale, color="#94a4b7", alpha=.07)
    for j, name in enumerate(MODELS):
        v = P@np.sqrt(total_force(name, grid, report))
        ax.plot(xx*scale, v, color=COLORS[name], lw=2.8 if name == "superposition" else 1.9,
                ls="--" if name in ("baryons", "nfw") else "-")
        y = .74-j*.08
        fig.text(.76, y, LABELS[name], color=COLORS[name], size=12)
        fig.text(.94, y-.027, f"{report['models'][name]['rms_kms']:.2f} km/s", ha="right", size=12, weight="bold")
    for used, color, marker in ((True, fg, "o"), (False, muted, "s")):
        select = mask == used
        ax.errorbar(obs["r_arcsec"][select]*scale, obs["v_los_kms"][select], yerr=obs["error_kms"][select],
                    fmt=marker, ms=4.3, color=color, mec=bg, capsize=3, lw=1, zorder=10)
    ax.set(xlabel="Projected major-axis radius [kpc]", ylabel="Observed line-of-sight rotation [km/s]",
           xlim=(0, 13.2), ylim=(0, 320))
    ax.grid(alpha=.12)
    fig.text(.76, .23, "RMS on the same six points", color=muted, size=10)
    fig.text(.08, .08, "Dutton et al. (2011), Table 4. Shaded inner region excluded, following the paper.", color=muted, size=10)
    fig.text(.08, .045, f"Fixed SPS baryons (mass offset {report['stellar_mass_offset_dex']:+.2f} dex); PSF + slit model. Signed terms: {report['models']['superposition']['parameter_count']}.", color=muted, size=10)
    fig.savefig(out/"j2141_rotation.png", dpi=180)
    fig.savefig(out/"j2141_rotation.pdf")
    plt.close(fig)

    fig, (left, right) = plt.subplots(1, 2, figsize=(13.5, 7.5), gridspec_kw={"width_ratios": [1.05, 1]})
    fig.subplots_adjust(left=.07, right=.96, top=.75, bottom=.23, wspace=.45)
    fig.text(.07, .93, "HORB  /  LENSING PREDICTION", color=COLORS["superposition"], size=13, weight="bold")
    fig.text(.07, .87, "SDSS J2141−0001", size=25, weight="bold")
    fig.text(.07, .82, "Mass parameters fixed by rotation; standard GR thin-lens projection", color=muted)
    maps = np.load(out/"j2141_lens_map.npz")
    dark = maps["kappa_dark"]
    limit = max(float(np.max(np.abs(dark))), .01)
    im = left.imshow(dark, origin="lower", extent=(-3.5, 3.5, -3.5, 3.5), cmap="RdBu_r",
                     norm=SymLogNorm(linthresh=max(limit*.01, 1e-4), vmin=-limit, vmax=limit))
    cr = curves["superposition"]
    phi = np.r_[cr["phi_radians"], 2*np.pi]
    radius = np.r_[cr["radius_kpc"], cr["radius_kpc"][0]]/scale
    if cr["closed_star_shaped_outer_curve"]:
        left.plot(radius*np.cos(phi), radius*np.sin(phi), color="#111111", lw=2.5)
    else:
        left.text(.5,.5,"Multiple critical branches\nSee lens-map data",transform=left.transAxes,
                  ha="center",va="center",color=fg,bbox={"facecolor":bg,"alpha":.85})
    left.set(xlabel="Major-axis offset [arcsec]", ylabel="Minor-axis offset [arcsec]",
             title="Projected HORB density + total critical curve")
    colorbar = fig.colorbar(im, ax=left, fraction=.046, pad=.025)
    colorbar.set_label(r"$\kappa_{\rm HORB}=\Sigma_{\rm HORB}/\Sigma_{\rm crit}$")
    obs_lens = lensing["source"]
    observed = obs_lens["circularized_einstein_radius_arcsec"]
    right.axvspan(observed-obs_lens["radius_error_minus_arcsec"], observed+obs_lens["radius_error_plus_arcsec"],
                  color=fg, alpha=.12)
    right.axvline(observed, color=fg, ls="--", lw=1.2)
    values = []
    alternative_values = [item["effective_einstein_radius_arcsec"] for item in lensing["alternative_rotation_fits"]
                          if item["effective_einstein_radius_arcsec"] is not None]
    if alternative_values:
        right.scatter(alternative_values, np.zeros(len(alternative_values)), color=COLORS["superposition"],
                      s=18, alpha=.25, zorder=3)
    for j, name in enumerate(MODELS):
        value = lensing["models"][name]["effective_einstein_radius_arcsec"]
        if value is not None:
            right.scatter(value, j, color=COLORS[name], s=60, zorder=5)
            right.annotate(f"{value:.3f}″", (value, j), xytext=(0, 10), textcoords="offset points",
                           ha="center", color=COLORS[name], size=10)
            values.append(value)
        else:
            right.text(.05,j,"multiple / missing critical branches",color=COLORS[name],size=9)
    right.set_yticks(range(len(MODELS)), [LABELS[n] for n in MODELS])
    right.set(xlabel="Area-equivalent Einstein radius [arcsec]", ylim=(5.6, -.7),
              xlim=(0, max(1.15, max(values)*1.12)), title="Frozen rotation-fit predictions")
    right.grid(axis="x", alpha=.12)
    fig.text(.07, .13, "White band: published SIE summary, 0.89 (+0.05 / −0.08) arcsec. Faint points: other near-best rotation fits; not a confidence interval.", color=muted, size=9.6)
    fig.text(.07, .09, "Conditional on fixed stellar masses, gas omission, orbital alignment and zero external shear/convergence.", color=muted, size=10)
    fig.text(.07, .05, "No observed-image pixel likelihood has been evaluated.", color=muted, size=10)
    fig.savefig(out/"j2141_lensing.png", dpi=180)
    fig.savefig(out/"j2141_lensing.pdf")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=ROOT/"reports/swells/j2141")
    p.add_argument("--max-components", type=int, default=4)
    p.add_argument("--beam-width", type=int, default=8)
    p.add_argument("--budget", type=int, default=1200)
    p.add_argument("--max-abs-weight", type=float, default=20.)
    p.add_argument("--stellar-mass-offset-dex", type=float, default=0.)
    p.add_argument("--nonnegative", action="store_true")
    p.add_argument("--exhaustive", action="store_true")
    p.add_argument("--fit-only", action="store_true")
    p.add_argument("--predict-only", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--validate-only", action="store_true")
    args = p.parse_args(argv)
    if (not 1 <= args.max_components <= 4 or args.beam_width < 1 or args.budget < 31
            or not np.isfinite(args.max_abs_weight) or not 0 < args.max_abs_weight <= 100
            or not np.isfinite(args.stellar_mass_offset_dex) or not -.5 <= args.stellar_mass_offset_dex <= .5
            or args.fit_only and args.predict_only):
        p.error("invalid fit controls or incompatible fit-only/predict-only options")
    meta, obs = read_inputs()
    if args.validate_only:
        print(f"Validated {len(obs)} J2141 measurements; {int(obs['use_for_fit'].sum())} fitted; five further SWELLS curves await numeric arrays.")
        return
    args.out_dir.mkdir(parents=True, exist_ok=True)
    path = args.out_dir/"rotation_fit.json"
    if args.predict_only:
        report = json.loads(path.read_text())
        if report["fit_inputs_sha256"] != fingerprints() or report["lensing_used_in_fit"]:
            raise ValueError("recorded fit has different inputs or used lensing; regenerate rotation fit")
    else:
        report = fit_rotation(meta, obs, args)
    rotation_outputs(report, meta, obs, args.out_dir)
    write_json(path, report)
    if args.fit_only:
        print(f"Wrote {path}")
        return
    lensing, curves = lens_predictions(report, meta, args.out_dir)
    if not args.no_plots:
        plot_results(report, lensing, curves, meta, obs, args.out_dir)
    print(f"Wrote rotation fits and held-out lensing predictions to {args.out_dir}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, RuntimeError) as exc:
        sys.exit(f"SWELLS: {exc}")
