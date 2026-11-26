"""Self-consistent density, equatorial force and observing geometry.

The baryons use the oblate chameleon profiles of SWELLS II, equations 3-8.
Mass normalizations are derived by integrating the densities, rather than
using the factor-of-two inconsistent total-mass expression in Appendix A8.
Hydrogenic forces retain all L=0,2,4 multipoles, exactly as the Rust backend.
"""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import (eval_legendre, gammainc, gammaincc, gamma,
                           genlaguerre, ndtr)
from horb_highlights import G, model_v2

C_KMS = 299792.458
ARCSEC_RAD = np.pi / (180 * 3600)
STATES = (("1s", 1, 0), ("2s", 2, 0), ("2p_z", 2, 1),
          ("3s", 3, 0), ("3p_z", 3, 1), ("3d_z2", 3, 2))


def distances(meta):
    zl, zs = meta["lens_redshift"], meta["source_redshift"]
    h0, om = meta["H0_kms_Mpc"], meta["Omega_m"]
    if not (0 < zl < zs and h0 > 0 and 0 < om < 1):
        raise ValueError("require 0 < lens z < source z, H0 > 0, 0 < Omega_m < 1")
    def chi(z):
        return C_KMS / h0 * quad(lambda x: 1 / np.sqrt(om*(1+x)**3 + 1-om), 0, z)[0] * 1000
    cl, cs = chi(zl), chi(zs)
    dl, ds, dls = cl/(1+zl), cs/(1+zs), (cs-cl)/(1+zs)
    return {"Dd_kpc": dl, "Ds_kpc": ds, "Dds_kpc": dls,
            "kpc_per_arcsec": dl*ARCSEC_RAD,
            "sigma_critical_msun_kpc2": C_KMS**2/(4*np.pi*G)*ds/(dl*dls)}


def baryon_components(meta, mass_offset_dex=0.0):
    distance = distances(meta)["kpc_per_arcsec"]
    ci = meta["cos_inclination"]
    result = []
    for c in meta["components"]:
        q = c["q_projected"]
        if not (0 < ci < q < 1 and c["core_arcsec"] > 0 and 0 < c["core_ratio"] < 1):
            raise ValueError("invalid oblate baryonic geometry")
        a = c["core_arcsec"] * distance
        result.append({"name": c["name"], "mass": 10**(c["log10_mass_msun"]+mass_offset_dex),
                       "q": q, "q3": np.sqrt((q*q-ci*ci)/(1-ci*ci)),
                       "a": a, "b": a/c["core_ratio"]})
    return result


def chameleon_v2(r, component):
    r = np.asarray(r, dtype=float)
    c = component
    e = np.sqrt(1-c["q3"]**2)
    f = e/np.arcsin(e) if e > 1e-6 else 1-e*e/6
    def shape(a):
        t = np.sqrt(r*r + e*e*a*a)
        # arctan(t/(q3*a))/t has a finite limit at the origin.
        ratio = np.divide(np.arctan(t/(c["q3"]*a)), t,
                          out=np.full_like(t, 1/(c["q3"]*a)), where=t != 0)
        return 1-f*a*ratio
    vnorm = 2*G*c["mass"]/(np.pi*f*(c["b"]-c["a"]))
    return np.maximum(vnorm*(shape(c["a"])-shape(c["b"])), 0)


def baryon_v2(r, components):
    return sum((chameleon_v2(r, c) for c in components), np.zeros_like(r, dtype=float))


@lru_cache(None)
def radial_polynomial(n, ell):
    if not (1 <= n <= 3 and 0 <= ell < n):
        raise ValueError("this analytic dictionary supports n <= 3")
    beta = 2/n
    p = genlaguerre(n-ell-1, 2*ell+1).coef[::-1]
    p = p * beta**np.arange(len(p))
    squared = np.polynomial.polynomial.polymul(p, p)
    norm = beta**3 * gamma(n-ell)/(2*n*gamma(n+ell+1)) * beta**(2*ell)
    return beta, norm*squared, np.arange(len(squared))+2*ell


def orbital_v2(r, n, ell, scale, mass=1e11):
    r = np.asarray(r, dtype=float)
    if scale <= 0 or np.any(r < 0) or not np.all(np.isfinite(r)):
        raise ValueError("finite nonnegative radius and positive scale required")
    x = np.maximum(r/scale, 1e-10)
    beta, coefficients, powers = radial_polynomial(n, ell)
    mu, weights = leggauss(8)
    angular_density = (2*ell+1)/(4*np.pi)*eval_legendre(ell, mu)**2
    result = np.zeros_like(x)
    for order in range(0, 2*ell+1, 2):
        angular = 2*np.pi*np.sum(weights*angular_density*eval_legendre(order, mu))*eval_legendre(order, 0)
        inside, outside = np.zeros_like(x), np.zeros_like(x)
        for coefficient, power in zip(coefficients, powers):
            k = power+order+3
            inside += coefficient*gamma(k)*gammainc(k, beta*x)/beta**k
            if order:
                k = power+2-order
                outside += coefficient*gamma(k)*gammaincc(k, beta*x)/beta**k
        result += angular*((order+1)*inside/x**(order+1) - order*x**order*outside)
    return np.where(r == 0, 0, G*mass/scale*result)


def dictionary(scales=(.5, 1., 2., 4., 8.), density_model="squared"):
    return [{"state": name, "n": n, "l": ell, "m": 0, "a0_star_kpc": a,
             "ref_mass_msun": 1e11, **({"density_model": "psi"} if density_model == "psi" else {})} for name, n, ell in STATES if density_model != "psi" or ell % 2 == 0 for a in scales]


def basis_matrix(r, terms):
    from horb_lobes import orbital_v2 as psi_v2
    return np.column_stack([(psi_v2 if t.get("density_model", "squared") == "psi" else orbital_v2)(
        r, t["n"], t["l"], t["a0_star_kpc"], t["ref_mass_msun"]) for t in terms]) if terms else np.zeros((len(r), 0))


def baseline_v2(model, r, params, reference=13.):
    r = np.asarray(r, dtype=float)
    return np.where(r == 0, 0, model_v2(model, np.maximum(r, 1e-8), params, reference))


def observing_matrix(radii_arcsec, grid_kpc, meta, step_arcsec=.04):
    """Linear map from intrinsic circular speeds to flux-weighted LOS speeds.

    Integrate an exponential emitting thin disk through a 1 arcsec slit and
    finite spatial bins, convolved with the Gaussian seeing PSF. The map is
    linear in velocity, NOT in v². Apply it after summing signed forces and
    taking the square root. Exact Gaussian rectangle integrals weight pixels.
    """
    rr = np.asarray(radii_arcsec, dtype=float)
    grid = np.asarray(grid_kpc, dtype=float)
    if (not np.all(np.isfinite(rr)) or np.any(rr < 0) or grid[0] != 0
            or np.any(np.diff(grid) <= 0) or not .005 <= step_arcsec <= .1):
        raise ValueError("invalid observing grid")
    sigma = meta["seeing_fwhm_arcsec"] / np.sqrt(8*np.log(2))
    half_bin, half_slit = meta["spatial_bin_arcsec"]/2, meta["slit_width_arcsec"]/2
    limit_x, limit_y = rr.max()+half_bin+6*sigma, half_slit+6*sigma
    # Symmetric midpoints preserve exactly zero central projected rotation.
    nx, ny = int(np.ceil(limit_x/step_arcsec)), int(np.ceil(limit_y/step_arcsec))
    x = (np.arange(2*nx)+.5-nx)*step_arcsec
    y = (np.arange(2*ny)+.5-ny)*step_arcsec
    xx, yy = np.meshgrid(x, y, indexing="ij")
    disk_r = np.hypot(xx, yy/meta["cos_inclination"])
    radius = disk_r*distances(meta)["kpc_per_arcsec"]
    if radius.max() > grid[-1]:
        raise ValueError("force grid too small for seeing/slit integration")
    upper = np.searchsorted(grid, radius, side="right")
    fraction = (radius-grid[upper-1])/(grid[upper]-grid[upper-1])
    emission = np.exp(-disk_r/meta["emission_scale_arcsec"])
    projection = np.sqrt(1-meta["cos_inclination"]**2)*xx/disk_r
    wy = ndtr((half_slit-y)/sigma)-ndtr((-half_slit-y)/sigma)
    result = np.zeros((len(rr), len(grid)))
    for j, center in enumerate(rr):
        wx = ndtr((center+half_bin-x)/sigma)-ndtr((center-half_bin-x)/sigma)
        flux = wx[:, None]*wy[None, :]*emission
        weighted = flux*projection/flux.sum()
        np.add.at(result[j], (upper-1).ravel(), (weighted*(1-fraction)).ravel())
        np.add.at(result[j], upper.ravel(), (weighted*fraction).ravel())
    return result
