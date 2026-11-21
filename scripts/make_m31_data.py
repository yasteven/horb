#!/usr/bin/env python3
"""Prepare the published M31 curve and a reproducible fixed baryonic baseline."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import i0e, i1e, k0e, k1e
from fit_signed_density import atomic_write, write_csv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/m31"
G = 4.3009e-6


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def observations(path):
    with Path(path).open(newline="") as f:
        rows = list(csv.DictReader(f))
    result = np.array([[float(row[k]) for k in ("R_kpc", "Vobs_kms", "eV_kms")] for row in rows])
    if (result.shape != (29, 3) or not np.all(np.isfinite(result))
            or not np.array_equal(result[:, 0], np.arange(8.5, 37))
            or np.any(result[:, 1:] <= 0)):
        raise ValueError("Corbelli Table 1 requires 29 positive measurements at 8.5,9.5,...,36.5 kpc")
    return result


def thin_disk_v2(r, mass, scale):
    """Geehan Eq. 7 / Freeman disk, not spherical enclosed-mass gravity."""
    r = np.asarray(r, dtype=float)
    if (not np.all(np.isfinite(r)) or np.any(r <= 0) or not np.isfinite(mass)
            or mass < 0 or not np.isfinite(scale) or scale <= 0):
        raise ValueError("disk requires positive finite radii/scale and nonnegative finite mass")
    y = r / (2*scale)
    # Opposite exponential scalings cancel in each I*K product.
    return 2*G*mass/scale*y*y*(i0e(y)*k0e(y)-i1e(y)*k1e(y))


def baryonic_components(r, model):
    p = model["parameters"]
    keys = {"bulge_mass_msun", "bulge_scale_kpc", "disk_mass_msun", "disk_scale_kpc", "black_hole_mass_msun"}
    if model["schema_version"] != 1 or set(p) != keys or any(not np.isfinite(p[k]) or p[k] <= 0 for k in keys):
        raise ValueError("invalid M31 baryonic-model parameters")
    r = np.asarray(r, dtype=float)
    if not np.all(np.isfinite(r)) or np.any(r <= 0):
        raise ValueError("baryons require positive finite radii")
    bulge = G*p["bulge_mass_msun"]*r/(r+p["bulge_scale_kpc"])**2
    disk = thin_disk_v2(r, p["disk_mass_msun"], p["disk_scale_kpc"])
    bh = G*p["black_hole_mass_msun"]/r
    return np.column_stack([bulge, disk, bh])


def prepare(out_dir, source=DATA/"corbelli2010_table1.csv", model_path=DATA/"geehan2006_model.json"):
    obs = observations(source)
    model = json.loads(Path(model_path).read_text())
    out_dir = Path(out_dir)
    # Dense force interpolation; exact observation radii are also present.
    r = np.unique(np.concatenate([np.arange(5, 2001, dtype=float)/20, obs[:, 0]]))
    components = baryonic_components(r, model)
    rc, baryons = out_dir/"rotation_curve.csv", out_dir/"baryons.csv"
    write_csv(rc, ["R_kpc", "Vobs_kms", "eV_kms"], [[f"{x:.1f}" for x in row] for row in obs])
    write_csv(baryons, ["r_kpc", "v2_bulge_kms2", "v2_disk_kms2", "v2_black_hole_kms2", "v2_baryons_kms2"],
              np.column_stack([r, components, components.sum(axis=1)]))
    metadata = {
        "schema_version": 1, "galaxy": "M31 (Andromeda)", "dataset_id": "m31_corbelli2010_geehan2006",
        "rotation_curve": {
            "reference": "Corbelli et al. (2010), A&A 511, A89, Table 1",
            "url": "https://arxiv.org/html/0912.4133v1", "doi": "10.1051/0004-6361/200913297",
            "distance_kpc": 785, "points": len(obs), "fit_range_kpc": [8.5, 36.5],
            "errors": "Published Eq. 5: half the difference between galaxy-side means plus their mean dispersion; not standard errors of the mean. No covariance matrix is supplied.",
            "geometry": "Published tilted-ring deprojection retained; no additional inclination correction",
            "excluded_inner_region": "Authors exclude R<8 kpc because of non-circular motions; no fabricated central points",
            "pdf_checked_sha256": "66aadd6a9d71570da56e95ac908caf37820c4ab187d1f329091c03d53ed08658"
        },
        "baryons": model,
        "inputs": {"observations": {"sha256": digest(rc)}, "baryons": {"sha256": digest(baryons)},
                   "source_table": {"name": Path(source).name, "sha256": digest(source)},
                   "model_parameters": {"name": Path(model_path).name, "sha256": digest(model_path)}},
        "notes": ["M31 is not in the SPARC 175-galaxy sample.",
                  "Fit errors use the tabulated error bars as diagonal weights; a small chi2/N is not a calibrated probability.",
                  "The deprojected curve includes a warp; HORB evaluates an axisymmetric flat-plane approximation.",
                  *model["notes"]]
    }
    atomic_write(out_dir/"dataset.json", lambda f: json.dump(metadata, f, indent=2, allow_nan=False))
    return metadata


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=ROOT/"reports/m31_inputs")
    args = p.parse_args()
    prepare(args.out_dir)
    print(f"M31: 29 Corbelli et al. (2010) observations, 8.5-36.5 kpc; published error bars retained")
    print("Fixed Geehan et al. (2006) bulge + thin disk + central mass; no separate gas component")
    print(f"wrote {args.out_dir}/rotation_curve.csv, baryons.csv, dataset.json")


if __name__ == "__main__":
    main()
