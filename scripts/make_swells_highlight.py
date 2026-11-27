#!/usr/bin/env python3
"""Rebuild the README's SWELLS figures from the recorded fits, including QI."""
import argparse
import csv
import json
from pathlib import Path
import shutil

from fit_swells import (ROOT, DATA, ROTATION_MODELS, read_inputs, rotation_outputs,
                        plot_rotation, write_json, ensure_psi_model, lens_predictions)
from compare_swells import load_family, compare, audit_and_map, plot_comparison


def lensing_density_report(comparison):
    """Use the same signed-psi candidate as the closest-lensing comparison dot."""
    selected = [row for row in comparison["rows"]
                if row["family"] == "psi" and row["selection"] == "closest_lensing"]
    if len(selected) != 1 or not selected[0]["selection_uses_lensing"]:
        raise ValueError("signed-psi density highlight requires one closest-lensing selection; "
                         "the rotation winner cannot substitute for it")
    return {**selected[0], "schema_version": 1, "galaxy": comparison["galaxy"],
            "density_model": "psi", "fit_inputs_sha256": comparison["fit_inputs_sha256"],
            "selection_statistic": comparison["lensing_selection_statistic"],
            "selection_scope": comparison["selection_scope"],
            "saved_candidates": comparison["pools"]["psi"]["saved_candidates"],
            "lensing_source": comparison["source"],
            "selection_source": "data/swells/reference/j2141_comparison/selection.json",
            "slice": {"plane": "intrinsic x-z", "y_kpc": 0.0}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=ROOT/"reports/swells/readme")
    args = parser.parse_args(argv)
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    signed = load_family(DATA/"reference/j2141", False)
    positive = load_family(DATA/"reference/j2141_nonnegative", True)
    report = signed[0]
    meta, obs = read_inputs()
    expected_obs = [{k: float(row[k]) for k in obs.dtype.names} for row in obs]
    if report["photometry"] != meta or report["observations"] != expected_obs:
        raise ValueError("recorded SWELLS geometry or observations differ from the source")
    if "psi" not in report["models"]:
        ensure_psi_model(report)
        lensing, _ = lens_predictions(report, meta, out)
        signed = (report, lensing)
        write_json(DATA/"reference/j2141/rotation_fit.json", report)
        write_json(DATA/"reference/j2141/lensing_predictions.json", lensing)
    rotation_outputs(report, meta, obs, out)
    plot_rotation(report, meta, obs, out)
    comparison = compare(signed, positive)
    lens_map = audit_and_map(comparison, signed, positive, out)
    density = lensing_density_report(comparison)
    # The intrinsic slice uses the lensing-selected psi weights, independently
    # of the rotation winner displayed in the separate rotation-curve panel.
    import numpy as np
    from horb_highlights import signed_density_xz
    from make_readme_highlights import GALAXIES, density_figure
    GALAXIES["j2141"] = {"title":"SDSS J2141−0001", "number":"03"}
    axis=np.linspace(-40,40,401)
    xx,zz=np.meshgrid(axis,axis)
    density_figure("j2141",axis,axis,signed_density_xz(xx,zz,density["components"]),
                   density,ROOT/"images/swells/j2141_psi_density_highlight.png")
    plot_comparison(comparison, lens_map, out)
    records = ROOT/"data/highlights"
    records.mkdir(parents=True, exist_ok=True)
    write_json(records/"j2141_psi_density_fit.json", density)
    # Persist the audited selections used for the figures, including psi.
    for directory in (out, DATA/"reference/j2141_comparison"):
        directory.mkdir(parents=True, exist_ok=True)
        write_json(directory/"selection.json", comparison)
        fields = ["family", "selection", "selection_uses_lensing", "rotation_rms_kms",
                  "rotation_chi2", "n_rotation", "effective_einstein_radius_arcsec"]
        with (directory/"summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            writer.writerows(comparison["rows"])
    scores = {name: report["numerical_audit"][name]["refined_scores"] for name in ROTATION_MODELS}
    scores["qi"] = {**report["models"]["qi"], **scores["qi"]}
    write_json(records/"j2141_comparison.json", {
        "galaxy": report["galaxy"], "fit_inputs_sha256": report["fit_inputs_sha256"],
        "fit_n": report["fit_n"], "observations": report["observations"],
        "observed_velocity_convention": report["observed_velocity_convention"],
        "photometry": meta, "baryons": report["baryons"], "scores": scores,
        "numerical_audit": report["numerical_audit"],
        "qi_lensing": "Not evaluated: the rotation prescription does not supply a lensing law.",
        "rotation_fit_source": "data/swells/reference/j2141/rotation_fit.json",
        "lensing_selection_source": "data/swells/reference/j2141_comparison/selection.json",
        "psi_density_fit_source": "data/highlights/j2141_psi_density_fit.json",
    })
    shutil.copyfile(out/"rotation_predictions.csv", records/"j2141_predictions.csv")
    images, pdfs = ROOT/"images/swells", ROOT/"plots"
    images.mkdir(parents=True, exist_ok=True)
    pdfs.mkdir(parents=True, exist_ok=True)
    for kind in ("rotation", "lensing"):
        shutil.copyfile(out/f"j2141_{kind}.png", images/f"j2141_{kind}.png")
        shutil.copyfile(out/f"j2141_{kind}.pdf", pdfs/f"j2141_{kind}_highlight.pdf")
    print(json.dumps({name: round(score["rms_kms"], 4) for name, score in scores.items()}))
    print("Updated J2141 rotation and lensing README figures from recorded fits.")


if __name__ == "__main__":
    main()
