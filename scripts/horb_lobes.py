"""Signed real psi mass sources, normalized by integral |psi| dV.

This is a phenomenological gravitational source, not a Born probability or
|sum(c psi)|². Global phase is fixed by positive Laguerre normalization and
the repository's real spherical harmonics. A sign flip is a weight flip.
Independent SciPy implementation audits the Rust multipole/density exports.
"""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_genlaguerre, eval_legendre, gamma, roots_genlaguerre

G = 4.3009e-6
LOBE_STATES = "1s,2s,3s,3d_z2"
MODEL_LABELS = {"squared": "HORB |ψ|²", "psi": "HORB signed ψ"}


def validate_state(n, ell):
    if not (1 <= n <= 10 and 0 <= ell < n and ell <= 2):
        raise ValueError("psi basis supports n<=10, l<=2")


def radial(r, n, ell):
    validate_state(n, ell)
    u = 2*np.asarray(r, dtype=float)/n
    norm = np.sqrt((2/n)**3*gamma(n-ell)/(2*n*gamma(n+ell+1)))
    return norm*np.exp(-u/2)*u**ell*eval_genlaguerre(n-ell-1, 2*ell+1, u)


@lru_cache(None)
def intervals(n, ell):
    validate_state(n, ell)
    roots = roots_genlaguerre(n-ell-1, 2*ell+1)[0]*n/2 if n-ell-1 else []
    return np.r_[0., roots, max(16*n*n, 120*n)]


@lru_cache(None)
def quadrature(order):
    return leggauss(order)


def moment(lo, hi, power, n, ell, absolute=False, order=128):
    lo, hi = np.broadcast_arrays(np.asarray(lo), np.asarray(hi))
    result = np.zeros_like(lo, dtype=float)
    nodes, weights = quadrature(order)
    for a, b in zip(intervals(n, ell)[:-1], intervals(n, ell)[1:]):
        left, right = np.maximum(lo, a), np.minimum(hi, b)
        half = np.maximum(right-left, 0)/2
        x = (left+half)[..., None]+half[..., None]*nodes
        # No contribution outside this interval; positive placeholder avoids
        # 0 * infinity for exterior moments with negative powers.
        x = np.where(half[..., None] > 0, x, 1.)
        values = radial(x, n, ell)
        if absolute:
            values = np.abs(values)
        result += half*np.sum(weights*values*x**power, axis=-1)
    return result


@lru_cache(None)
def normalization(n, ell):
    angular_abs = {0: np.sqrt(4*np.pi), 1: np.sqrt(3*np.pi),
                   2: 4*np.sqrt(5*np.pi)/(3*np.sqrt(3))}[ell]
    return float(angular_abs*moment(0, intervals(n, ell)[-1], 2, n, ell, True))


def net_mass_fraction(n, ell):
    return float(np.sqrt(4*np.pi)*moment(0, intervals(n, ell)[-1], 2, n, ell)/normalization(n, ell)) if ell == 0 else 0.


def density(r, mu, n, ell, scale, mass):
    return (mass/scale**3/normalization(n, ell)*radial(np.asarray(r)/scale, n, ell)
            *np.sqrt((2*ell+1)/(4*np.pi))*eval_legendre(ell, mu))


def orbital_v2(r, n, ell, scale, mass=1e11):
    r = np.asarray(r, dtype=float)
    validate_state(n, ell)
    if scale <= 0 or mass < 0 or not np.isfinite([scale, mass]).all() or not np.isfinite(r).all() or np.any(r < 0):
        raise ValueError("invalid signed-lobe scale, mass or radius")
    if ell % 2:
        return np.zeros_like(r)
    x = np.maximum(r/scale, 1e-10)
    inside = moment(0, x, ell+2, n, ell)
    outside = moment(x, intervals(n, ell)[-1], 1-ell, n, ell) if ell else 0.
    angular = np.sqrt(4*np.pi/(2*ell+1))*eval_legendre(ell, 0)/normalization(n, ell)
    force = G*mass/scale*angular*((ell+1)*inside/x**(ell+1)-ell*x**ell*outside)
    return np.where(r == 0, 0, force)


def mass_fields(term):
    """Keep amplitude mass distinct from net mass for zero-monopole modes.

    Positive/negative masses here sum the lobes of individual components;
    after overlap, the positive/negative masses of the total field can differ.
    """
    amplitude = float(term["weight"]*term["ref_mass_msun"])
    fraction = net_mass_fraction(term["n"], term["l"]) if term.get("density_model", "squared") == "psi" else 1.
    net = amplitude*fraction
    return {"amplitude_mass_msun": amplitude, "absolute_component_mass_msun": abs(amplitude),
            "net_mass_fraction": fraction, "signed_mass_msun": net,
            "positive_lobe_mass_msun": (abs(amplitude)+net)/2,
            "negative_lobe_mass_msun": (net-abs(amplitude))/2}
