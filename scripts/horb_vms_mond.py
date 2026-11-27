"""VMS's downloadable curve prescription and the simple algebraic MOND law.

VMS is implemented as the author's script profile, not as a claim that the
PDF's distinct density, potential and calibration expressions are equivalent.
"""
import numpy as np
from scipy.optimize import brentq

from horb_qi import KPC_M

VMS_URL = "https://www.vms-institute.org/dark-matter/VMS_DARK_Matter_script.py"
VMS_SHA256 = "80d6b5301b17f89e5a56daf7458a8ce97e13b8768b95d1f0aa46a6271976e74f"
MOND_URL = "https://doi.org/10.12942/lrr-2012-10"
MOND_A0_M_S2 = 1.2e-10


def vms_v2(r_kpc, amplitude_kms2, scale_kpc):
    """Geometric-halo v² = A*r/(r+s), where A = v_infinity², not a density."""
    r = np.asarray(r_kpc, dtype=float)
    if (not np.isfinite(r).all() or np.any(r < 0)
            or not np.isfinite([amplitude_kms2, scale_kpc]).all()
            or amplitude_kms2 < 0 or scale_kpc <= 0):
        raise ValueError("VMS requires nonnegative radii/amplitude and positive scale")
    return amplitude_kms2 * r/(r+scale_kpc)


def fit_vms(r, observed, errors, baryons_v2):
    """Fit one amplitude in velocity-space chi²; s = max(observed r)/3.

    For fixed s, chi² is convex in A on the nonnegative total-force domain.
    Bracket its monotone derivative rather than depending on an initial guess.
    Preserve signed baryons and original errors instead of the source script's
    unweighted stellar components, squaring of gas signs and 0.5 km/s floor.
    """
    r, v, e, b = [np.asarray(a, dtype=float) for a in (r, observed, errors, baryons_v2)]
    if (r.ndim != 1 or not len(r) or any(a.shape != r.shape for a in (v, e, b))
            or not np.isfinite([r, v, e, b]).all() or np.any(r <= 0)
            or np.any(v < 0) or np.any(e <= 0)):
        raise ValueError("invalid VMS fit data")
    scale = float(r.max()/3)
    shape = r/(r+scale)
    lower = max(0., float(np.max(-b/shape)))

    def derivative(amplitude):
        total = np.maximum(b+amplitude*shape, 0.)
        ratio = np.divide(v, np.sqrt(total), out=np.where(v == 0, 0., np.inf), where=total > 0)
        return float(np.sum(shape*(1-ratio)/e**2))

    if derivative(lower) >= 0:
        amplitude = lower
    else:
        upper = max(1., 2*lower, float(np.max(v*v)/shape.min()))
        while derivative(upper) < 0:
            upper *= 2
            if not np.isfinite(upper):
                raise ValueError("could not bracket finite VMS amplitude")
        amplitude = brentq(derivative, lower, upper, xtol=1e-8)
    prediction = np.sqrt(np.maximum(b+amplitude*shape, 0.))
    delta = prediction-v
    return {"rms_kms": float(np.sqrt(np.mean(delta**2))), "chi2": float(np.sum((delta/e)**2)),
            "n": len(r), "source": VMS_URL, "source_sha256": VMS_SHA256,
            "implementation": "VMS downloadable script's rotation profile; HORB fixed baryons and errors",
            "parameter_count": 1, "fitted": True,
            "amplitude_kms2": float(amplitude), "v_infinity_kms": float(np.sqrt(amplitude)),
            "scale_kpc": scale, "scale_rule": "max(observed_radius_kpc)/3; fixed, not fitted",
            "amplitude_lower_bound_kms2": lower,
            "bound_reached": bool(np.isclose(amplitude, lower, rtol=1e-10, atol=1e-8)),
            "optimizer": "Convex velocity-space chi2: bracketed amplitude derivative root or lower boundary"}


def mond_velocity(r_kpc, baryons_v2_kms2):
    """Simple algebraic MOND: g = (g_bar + sqrt(g_bar²+4*a0*g_bar))/2.

    Fixed a0; no external-field effect or modified-Poisson field solver.
    Unlike the VMS download's deep-MOND-only control, this has both limits.
    """
    r, b = np.broadcast_arrays(np.asarray(r_kpc, dtype=float),
                               np.asarray(baryons_v2_kms2, dtype=float))
    if not np.isfinite(r).all() or np.any(r <= 0):
        raise ValueError("MOND requires finite positive radii")
    if not np.isfinite(b).all() or np.any(b < 0):
        raise ValueError("algebraic MOND requires finite nonnegative total baryonic v^2")
    a0r = MOND_A0_M_S2*KPC_M/1e6*r
    return np.sqrt(.5*(b+np.sqrt(b)*np.sqrt(b+4*a0r)))


def mond_metadata():
    return {"source": MOND_URL, "equations": [42, 43], "interpolation": "simple algebraic",
            "a0_m_s2": MOND_A0_M_S2, "parameter_count": 0, "fitted": False,
            "external_field": "omitted", "geometry": "algebraic rotation-curve approximation; no field solver"}
