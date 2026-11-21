"""Shared baryons, Einasto labels and scores for fixed halo comparisons."""

import sys
import numpy as np

G_KPC = 4.3009e-6


def baryon_velocity_component(kind, mass, scale, r):
    # Existing HORB baryonic prescription, unchanged by the Einasto addition.
    r = np.asarray(r)
    if kind in ("disk", "exponential_disk", "exponential-disk"):
        x = r / scale
        menc = mass * (1.0 - np.exp(-x) * (1.0 + x))
    elif kind in ("bulge", "hernquist", "hernquist_bulge", "hernquist-bulge"):
        menc = mass * r**2 / (r + scale) ** 2
    else:
        raise ValueError(f"unknown baryon component kind: {kind}")
    return np.sqrt(np.where(r > 0, G_KPC * menc / r, 0.0))


def load_baryon_velocity_squared(path, r):
    rows = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None))
    v2 = np.zeros_like(r, dtype=float)
    for row in rows:
        v = baryon_velocity_component(str(row["kind"]), float(row["mass_msun"]), float(row["scale_kpc"]), r)
        v2 += v**2
    return v2


def has_einasto(data):
    if "v_einasto_kms" in data.dtype.names:
        return True
    print("This archived CSV has no Einasto curve; regenerate the comparison to include it.", file=sys.stderr)
    return False


def einasto_label(data):
    fields = ("einasto_alpha", "einasto_r_minus2_kpc")
    if not all(name in data.dtype.names for name in fields):
        return "Einasto"
    alpha, scale = (float(np.atleast_1d(data[name])[0]) for name in fields)
    return rf"Einasto ($\alpha={alpha:g}$, $r_{{-2}}={scale:g}$ kpc)"


def score(name, r_model, v_model, observed, r_min, r_max):
    """Score the common radius interval; never clamp observations past a model."""
    r_model, v_model = np.asarray(r_model), np.asarray(v_model)
    if (r_model.size < 2 or r_model.shape != v_model.shape
            or not np.all(np.isfinite(r_model)) or np.any(np.diff(r_model) <= 0)
            or not np.all(np.isfinite(v_model))):
        raise ValueError("model requires finite velocities and at least two strictly increasing radii")
    if not (np.isfinite(r_min) and np.isfinite(r_max) and 0 <= r_min < r_max):
        raise ValueError("score range must be finite with 0 <= r_min < r_max")
    lower, upper = max(r_min, r_model[0]), min(r_max, r_model[-1])
    radii = np.atleast_1d(observed["R_kpc"])
    if not np.all(np.isfinite(radii)):
        raise ValueError("observed radii must be finite")
    mask = (radii >= lower) & (radii <= upper)
    if not np.any(mask):
        raise ValueError("no observations in the common model/score radius range")
    v_col, e_col = ("Vobs_kms", "eV_kms") if "Vobs_kms" in observed.dtype.names else ("v_kms", "v_err_kms")
    v_obs, err = (np.atleast_1d(observed[col])[mask] for col in (v_col, e_col))
    if not np.all(np.isfinite(v_obs)) or not np.all(np.isfinite(err) & (err > 0)):
        raise ValueError("selected observations require finite velocities and positive finite errors")
    residual = np.interp(radii[mask], r_model, v_model) - v_obs
    n = residual.size
    rms = np.sqrt(np.mean(residual**2))
    chi2 = np.sum((residual / err)**2)
    return name, lower, upper, n, rms, chi2, chi2 / n
