"""Offline, checked SPARC inputs. Preserve published velocities and force signs."""
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from fit_signed_density import atomic_write, write_csv

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/sparc"
G = 4.3009e-6  # kpc (km/s)^2 / Msun, same as orbital_basis
STATES = ("1s", "2s", "2p_z", "3s", "3p_z", "3d_z2")
CATALOG_FIELDS = (
    "galaxy", "hubble_type", "distance_mpc", "distance_error_mpc", "distance_method",
    "inclination_deg", "inclination_error_deg", "luminosity_1e9_lsun",
    "luminosity_error_1e9_lsun", "effective_radius_kpc", "effective_surface_brightness",
    "disk_scale_kpc", "disk_central_surface_brightness", "hi_mass_1e9_msun",
    "hi_radius_kpc", "v_flat_kms", "v_flat_error_kms", "quality", "references",
)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def save_json(path, data):
    atomic_write(path, lambda f: json.dump(data, f, indent=2, allow_nan=False))


def table_body(path):
    lines = Path(path).read_text().splitlines()
    separators = [i for i, line in enumerate(lines) if len(line) > 20 and set(line) == {"-"}]
    if not separators:
        raise ValueError(f"{path}: missing MRT header terminator")
    return [line.split() for line in lines[separators[-1]+1:] if line.strip()]


def load_data(source=DATA):
    source = Path(source)
    manifest = json.loads((source/"sources.json").read_text())
    expected = {"SPARC_Lelli2016c.mrt", "MassModels_Lelli2016c.mrt"}
    if manifest.get("schema_version") != 1 or set(manifest.get("files", {})) != expected:
        raise ValueError("unsupported SPARC source manifest")
    for name, info in manifest["files"].items():
        if digest(source/name) != info["sha256"]:
            raise ValueError(f"SPARC source checksum mismatch: {name}")
    catalog = {}
    # Actual public table has spaces between columns beyond its stated fixed
    # widths; all 19 fields are unambiguous whitespace-separated tokens.
    for fields in table_body(source/"SPARC_Lelli2016c.mrt"):
        if len(fields) != len(CATALOG_FIELDS):
            raise ValueError("malformed SPARC catalog row")
        name = fields[0]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in catalog:
            raise ValueError(f"invalid/duplicate galaxy: {name}")
        values = [name, *map(float, fields[1:-1]), fields[-1]]
        row = dict(zip(CATALOG_FIELDS, values))
        for key in ("hubble_type", "distance_method", "quality"):
            if not row[key].is_integer():
                raise ValueError(f"{name}: noninteger {key}")
            row[key] = int(row[key])
        if (not np.isfinite(values[1:-1]).all() or row["distance_mpc"] <= 0
                or row["disk_scale_kpc"] <= 0 or row["quality"] not in (1, 2, 3)
                or not 0 < row["inclination_deg"] <= 90):
            raise ValueError(f"{name}: invalid catalog metadata")
        catalog[name] = row
    curves = {name: [] for name in catalog}
    for fields in table_body(source/"MassModels_Lelli2016c.mrt"):
        if len(fields) != 10 or fields[0] not in catalog:
            raise ValueError("malformed/unknown SPARC mass-model row")
        name = fields[0]
        distance, *values = map(float, fields[1:])
        if distance != catalog[name]["distance_mpc"]:
            raise ValueError(f"{name}: catalog and rotation-curve distances disagree")
        curves[name].append(values)
    for name, rows in curves.items():
        a = np.asarray(rows)
        if (a.ndim != 2 or a.shape[1] != 8 or len(a) < 3 or not np.isfinite(a).all()
                or np.any(a[:, 0] <= 0) or np.any(np.diff(a[:, 0]) <= 0)
                or np.any(a[:, 1] < 0) or np.any(a[:, 2] <= 0)):
            raise ValueError(f"{name}: invalid rotation curve")
        curves[name] = a
    if (len(catalog) != 175 or manifest["galaxies"] != 175
            or sum(map(len, curves.values())) != 3391 or manifest["rotation_curve_points"] != 3391):
        raise ValueError("SPARC snapshot must contain all 175 galaxies and 3391 points")
    return catalog, curves, manifest


def select_galaxies(value, catalog):
    names = sorted(catalog) if value == "all" else value.split(",")
    if not names or len(names) != len(set(names)) or any(name not in catalog for name in names):
        raise ValueError("--galaxies requires 'all' or distinct exact SPARC names separated by commas")
    return sorted(names)


def mass_to_light(disk, bulge):
    if not np.isfinite([disk, bulge]).all() or min(disk, bulge) < 0:
        raise ValueError("stellar mass-to-light ratios must be finite and nonnegative")


def baryonic_components(curve, ml_disk, ml_bulge):
    mass_to_light(ml_disk, ml_bulge)
    velocities = curve[:, 3:6]
    return velocities * np.abs(velocities) * [1.0, ml_disk, ml_bulge]


def prepare_galaxy(out_dir, catalog_row, curve, manifest, ml_disk=0.5, ml_bulge=0.7):
    out_dir = Path(out_dir)
    components = baryonic_components(curve, ml_disk, ml_bulge)
    rc, baryons = out_dir/"rotation_curve.csv", out_dir/"baryons.csv"
    write_csv(rc, ["R_kpc", "Vobs_kms", "eV_kms"], curve[:, :3])
    write_csv(baryons, ["r_kpc", "v2_gas_kms2", "v2_disk_kms2", "v2_bulge_kms2", "v2_baryons_kms2"],
              np.column_stack([curve[:, 0], components, components.sum(axis=1)]))
    info = {
        "schema_version": 1, "galaxy": catalog_row["galaxy"], "dataset_id": "sparc_lelli2016",
        "catalog": catalog_row, "source": manifest,
        "rotation_curve": {
            "reference": manifest["reference"], "doi": manifest["doi"],
            "original_reference_keys": catalog_row["references"],
            "points": len(curve), "fit_range_kpc": [float(curve[0, 0]), float(curve[-1, 0])],
            "distance_mpc": catalog_row["distance_mpc"],
            "errors": "Published random errors; inclination/distance systematics are not added or fitted. Diagonal weighting; no covariance supplied.",
            "geometry": "Published deprojection retained; no additional inclination correction",
        },
        "baryons": {
            "ml_disk": ml_disk, "ml_bulge": ml_bulge, "ml_units": "Msun/Lsun at 3.6 microns",
            "plot_label": f"SPARC gas + disk (M/L={ml_disk:g}) + bulge (M/L={ml_bulge:g})",
            "formula": "Vgas*abs(Vgas) + ml_disk*Vdisk*abs(Vdisk) + ml_bulge*Vbul*abs(Vbul)",
            "helium": "Already included in published Vgas; no additional 1.33 multiplier",
            "interpolation": "Linear in signed v^2 within measured radii only; no extrapolation",
        },
        "inputs": {"observations": {"sha256": digest(rc)}, "baryons": {"sha256": digest(baryons)}},
        "notes": ["All published points retained, including quality 3 and low-inclination galaxies.",
                  "Negative component velocity denotes outward radial force, not negative gas mass.",
                  "Fit quality and inferred signed densities depend on these fixed baryons and the chosen basis.",
                  "No jet classifications or correlations are included."],
    }
    save_json(out_dir/"dataset.json", info)
    return info


def inventory(out_dir, catalog, curves, names):
    write_csv(Path(out_dir)/"catalog.csv", [*CATALOG_FIELDS, "points", "r_min_kpc", "r_max_kpc"],
              [[*[catalog[name][k] for k in CATALOG_FIELDS], len(curves[name]),
                curves[name][0, 0], curves[name][-1, 0]] for name in names])
