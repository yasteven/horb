"""Common-data fits for the README figures; equations match orbital_basis.

The 3d_z2 force is the exact L=0,2,4 multipole integral, evaluated with
incomplete gamma functions. It is NOT a spherical enclosed-mass approximation.
Classical halo masses are parameterized by M(<r_ref), not total mass.
"""
import numpy as np
from scipy.optimize import differential_evolution, least_squares
from scipy.special import gammainc, gammaincc, gamma, eval_genlaguerre, eval_legendre

G = 4.3009e-6
MODELS = ("superposition", "single", "nfw", "burkert", "einasto", "baryons")
LABELS = {"superposition": "HORB superposition", "single": "HORB single orbital",
          "nfw": "NFW", "burkert": "Burkert", "einasto": "Einasto", "baryons": "Baryons alone"}


def dz2_v2(r, scale, mass):
    """Unsoftened disk-plane v² for M |R_32 Y_20|², including exterior mass."""
    x = np.asarray(r, dtype=float) / scale
    if np.any(x <= 0) or scale <= 0 or mass < 0:
        raise ValueError("positive radii/scale and nonnegative mass required")
    beta, c = 2 / 3, 8 / 98415
    result = np.zeros_like(x)
    for ell, angular in [(0, 1), (2, -1 / 7), (4, 3 / 28)]:
        inside = c * gamma(ell + 7) * gammainc(ell + 7, beta*x) / beta**(ell + 7)
        outside = c * gamma(6 - ell) * gammaincc(6 - ell, beta*x) / beta**(6 - ell)
        result += angular * ((ell+1)*inside/x**(ell+1) - ell*x**ell*outside)
    return G * mass / scale * result


def halo_mass_shape(model, x, alpha=0.18):
    x = np.asarray(x, dtype=float)
    if model == "nfw":
        # Stable at small radii: x²/2 - 2x³/3 + 3x⁴/4 - ...
        series = sum((-1)**k * (k-1)/k * x**k for k in range(2, 9))
        return np.where(x < 1e-3, series, np.log1p(x) - x/(1+x))
    if model == "burkert":
        series = 4*x**3/3 - x**4 + 4*x**7/7 - x**8/2 + 4*x**11/11
        return np.where(x < 1e-3, series, 2*np.log1p(x) + np.log1p(x*x) - 2*np.arctan(x))
    if model == "einasto":
        return gammainc(3/alpha, 2/alpha * x**alpha)
    raise ValueError(f"unknown halo: {model}")


def model_v2(model, r, parameters, reference_radius):
    mass, scale = 10**parameters[0], 10**parameters[1]
    if model == "single":
        return dz2_v2(r, scale, mass)
    alpha = 10**parameters[2] if model == "einasto" else 0.18
    ratio = halo_mass_shape(model, np.asarray(r)/scale, alpha) / halo_mass_shape(model, reference_radius/scale, alpha)
    return G * mass / np.asarray(r) * ratio


def metrics(predicted, observed, errors):
    delta = np.asarray(predicted) - observed
    return {"rms_kms": float(np.sqrt(np.mean(delta**2))),
            "chi2": float(np.sum((delta/errors)**2)), "n": len(delta)}


def fit_baseline(model, r, observed, errors, baryons_v2, reference_radius):
    """Deterministic bounded global search followed by velocity-space LS."""
    bounds = [(8., 13.), (np.log10(.05), np.log10(50 if model == "single" else 500))]
    if model == "einasto":
        bounds.append((np.log10(.03), np.log10(2.)))
    def residual(p):
        total = baryons_v2 + model_v2(model, r, p, reference_radius)
        if np.any(total <= 0) or not np.all(np.isfinite(total)):
            raise ValueError("invalid total circular force in baseline fit")
        return (np.sqrt(total)-observed)/errors
    global_fit = differential_evolution(lambda p: np.sum(residual(p)**2), bounds,
                                        seed=42, popsize=12, maxiter=600, tol=1e-10, polish=False)
    local = least_squares(residual, global_fit.x, bounds=np.array(bounds).T,
                          jac="3-point", ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=3000)
    if not local.success or not np.all(np.isfinite(local.fun)):
        raise ValueError(f"{model} fit did not converge: {local.message}")
    params = local.x
    at_bound = any(min(p-lo, hi-p) < 1e-5 for p, (lo, hi) in zip(params, bounds))
    prediction = np.sqrt(baryons_v2 + model_v2(model, r, params, reference_radius))
    return {**metrics(prediction, observed, errors), "log10_parameters": params.tolist(),
            "parameter_names": ["total_mass_msun" if model == "single" else "mass_within_reference_msun", "scale_kpc"]
                               + (["alpha"] if model == "einasto" else []),
            "parameters": (10**params).tolist(), "log10_bounds": bounds,
            "reference_radius_kpc": reference_radius, "parameter_count": len(params),
            "bound_reached": bool(at_bound), "optimizer": "differential_evolution(seed=42) + least_squares",
            "global_converged": bool(global_fit.success), "local_converged": bool(local.success)}


def signed_density_xz(x, z, components):
    """Independent density evaluation to audit the Rust CSV before plotting."""
    radius = np.hypot(x, z)
    mu = np.divide(z, radius, out=np.ones_like(radius), where=radius != 0)
    rho = np.zeros_like(radius)
    for term in components:
        n, ell, m, a = (term[k] for k in ("n", "l", "m", "a0_star_kpc"))
        if m != 0:
            raise ValueError("README dictionary requires axisymmetric m=0 orbitals")
        u = 2*radius/(n*a)
        norm2 = (2/(n*a))**3 * gamma(n-ell)/(2*n*gamma(n+ell+1))
        radial2 = norm2*np.exp(-u)*u**(2*ell)*eval_genlaguerre(n-ell-1, 2*ell+1, u)**2
        angular2 = (2*ell+1)/(4*np.pi)*eval_legendre(ell, mu)**2
        rho += term["weight"]*term["ref_mass_msun"]*radial2*angular2
    return rho


def load_density(path, components):
    data = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))
    x, z = np.unique(data["x_kpc"]), np.unique(data["z_kpc"])
    if len(data) != len(x)*len(z) or min(len(x), len(z)) < 3:
        raise ValueError("density must be a complete rectangular grid")
    rows = data[np.lexsort((data["x_kpc"], data["z_kpc"]))]
    xx, zz = np.meshgrid(x, z)
    if not np.array_equal(rows["x_kpc"], xx.ravel()) or not np.array_equal(rows["z_kpc"], zz.ravel()):
        raise ValueError("density grid has missing or duplicate coordinates")
    rho = rows["rho_signed_msun_per_kpc3"].reshape(xx.shape)
    expected = signed_density_xz(xx, zz, components)
    if not np.all(np.isfinite(rho)) or not np.allclose(rho, expected, rtol=2e-10, atol=max(np.max(np.abs(expected))*1e-12, 1e-12)):
        raise ValueError("density CSV does not match the selected superposition weights")
    return x, z, rho
