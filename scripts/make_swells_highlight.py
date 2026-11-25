#!/usr/bin/env python3
"""Rebuild the README's SWELLS figures from the recorded fits, including QI."""
import argparse
import json
from pathlib import Path
import shutil

from fit_swells import (ROOT, DATA, ROTATION_MODELS, read_inputs, rotation_outputs,
                        plot_rotation, write_json)
from compare_swells import load_family, compare, audit_and_map, plot_comparison


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
    rotation_outputs(report, meta, obs, out)
    plot_rotation(report, meta, obs, out)
    comparison = compare(signed, positive)
    lens_map = audit_and_map(comparison, signed, positive, out)
    plot_comparison(comparison, lens_map, out)
    records = ROOT/"data/highlights"
    records.mkdir(parents=True, exist_ok=True)
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
