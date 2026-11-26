"""Thin-lens predictions under the weak-field GR mass-to-light-deflection law.

Signed coefficients are preserved. This assumes the fitted effective density
sources both metric potentials equally; it does not derive a relativistic
completion of HORB. The background source redshift fixes Sigma_crit.

Axisymmetric m=0 hydrogenic densities with l<=2 project into exactly the
cos(0 phi), cos(2 phi), cos(4 phi) modes. Their 2D Poisson solution can therefore
be integrated radially without pixel softening, sphericalizing the density,
or imposing a finite FFT box. Chameleon baryons have analytic deflections.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import eval_genlaguerre, gamma, gammaln
from horb_highlights import halo_mass_shape
from .quadrature import cumulative_simpson_1d


def orbital_surface_modes(radii, term, cos_i, order=192):
    r = np.asarray(radii, dtype=float)[:, None]
    n, ell, a, mass = (term[k] for k in ("n", "l", "a0_star_kpc", "ref_mass_msun"))
    nodes, weights = leggauss(order)
    # s = r*sinh(t) makes the integral resolve the tangent point and decades
    # of line-of-sight distance, including very compact density components.
    t, w = (nodes+1)*20, weights*20
    s, radius = r*np.sinh(t), r*np.cosh(t)
    u = 2*radius/(n*a)
    radial = (2/(n*a))**3*gamma(n-ell)/(2*n*gamma(n+ell+1))
    radial = radial*np.exp(-u)*u**(2*ell)*eval_genlaguerre(n-ell-1, 2*ell+1, u)**2
    lobes = term.get("density_model", "squared") == "psi"
    if lobes:
        if ell % 2 or term.get("m", 0) != 0:
            raise ValueError("psi lens projection supports aligned even-l, m=0 modes")
        from horb_lobes import radial as psi_radial, normalization
        radial = mass/a**3*psi_radial(radius/a, n, ell)/normalization(n, ell)*np.sqrt((2*ell+1)/(4*np.pi))
    else:
        radial = mass*radial*(2*ell+1)/(4*np.pi)
    si2, ci2 = 1-cos_i*cos_i, cos_i*cos_i
    surfaces = []
    for sin_phi2 in (0., .5, 1.):
        y2 = r*r*sin_phi2
        mu2 = (y2*si2+s*s*ci2)/radius**2
        if ell == 0:
            angular = 1.
        elif lobes:
            angular = (3*mu2-1)/2
        elif ell == 1:
            angular = mu2
        else:
            mu4 = (y2*y2*si2**2+6*y2*si2*s*s*ci2+s**4*ci2**2)/radius**4
            angular = (9*mu4-6*mu2+1)/4
        surfaces.append(2*np.sum(radial*angular*radius*w, axis=1))
    x, diagonal, y = surfaces
    return np.array([(x+y+2*diagonal)/4, (x-y)/2, (x+y-2*diagonal)/4])


def spherical_surface(radii, model, params, reference=13., order=192):
    mass, scale = 10**params[0], 10**params[1]
    alpha = 10**params[2] if model == "einasto" else .18
    nodes, weights = leggauss(order)
    t, w = (nodes+1)*20, weights*20
    radius = np.asarray(radii)[:, None]*np.cosh(t)
    x = radius/scale
    shape = halo_mass_shape(model, reference/scale, alpha)
    if model == "nfw":
        rho = mass/(4*np.pi*scale**3*shape) / (x*(1+x)**2)
    elif model == "burkert":
        rho = mass/(np.pi*scale**3*shape) / ((1+x)*(1+x*x))
    elif model == "einasto":
        lognorm = (np.log(mass)-np.log(4*np.pi)-3*np.log(scale)+np.log(alpha)
                   -3/alpha*np.log(alpha/2)-gammaln(3/alpha)-np.log(shape))
        rho = np.exp(lognorm-2/alpha*x**alpha)
    else:
        raise ValueError("spherical lens requires NFW, Burkert, or Einasto")
    return 2*np.sum(rho*radius*w, axis=1)


class FourierLens:
    def __init__(self, radii, kappa_modes):
        r = np.asarray(radii)
        kappa = np.asarray(kappa_modes)
        if (kappa.shape != (3, len(r)) or np.any(r <= 0) or np.any(np.diff(r) <= 0)
                or not np.all(np.isfinite(kappa))):
            raise ValueError("invalid projected density modes")
        self.rmin, self.rmax = r[0], r[-1]
        t = np.log(r)
        self.coefficients = []
        for index, m in enumerate((0, 2, 4)):
            k = kappa[index]
            inside = cumulative_simpson_1d(k*r**(m+2), x=t, initial=0)
            # Central monopole mass is finite; retain its small integral.
            if m == 0:
                inside += k[0]*r[0]**2/2
                a = 2*inside/r
                b = np.zeros_like(r)
                derivative = 2*k-a/r
            else:
                outside = cumulative_simpson_1d((k*r**(2-m))[::-1], x=-t[::-1], initial=0)[::-1]
                inner, outer = inside/r**(m+1), outside*r**(m-1)
                a, b = inner-outer, inner+outer
                derivative = 2*k-(m+1)*inner/r-(m-1)*outer/r
            self.coefficients.append(CubicSpline(t, np.array([a, b, derivative]), axis=1))

    def field(self, x, y):
        """Return alpha_x, alpha_y [kpc], Hxx, Hxy, Hyy [dimensionless]."""
        x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
        radius = np.hypot(x, y)
        if np.any(radius > self.rmax):
            raise ValueError("requested lens position exceeds the radial integration grid")
        r = np.maximum(radius, self.rmin*2)
        phi = np.arctan2(y, x)
        ar, ap, hrr, hrp, hpp = [np.zeros_like(r) for _ in range(5)]
        for m, spline in zip((0, 2, 4), self.coefficients):
            a, b, derivative = spline(np.log(r))
            co, si = np.cos(m*phi), np.sin(m*phi)
            ar += a*co
            ap += b*si
            hrr += derivative*co
            hrp -= (m*a+b)/r*si
            hpp += (a+m*b)/r*co
        co, si = np.cos(phi), np.sin(phi)
        ax, ay = ar*co-ap*si, ar*si+ap*co
        hxx = co*co*hrr-2*co*si*hrp+si*si*hpp
        hyy = si*si*hrr+2*co*si*hrp+co*co*hpp
        hxy = co*si*(hrr-hpp)+(co*co-si*si)*hrp
        return np.array([np.where(radius == 0, 0, ax), np.where(radius == 0, 0, ay), hxx, hxy, hyy])


def chameleon_field(x, y, c, sigma_critical):
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    q = c["q"]
    e = np.sqrt(1-q*q)
    strength = c["mass"]/(np.pi*sigma_critical*(c["b"]-c["a"]))
    def nie(a):
        psi = np.sqrt(q*q*(a*a+x*x)+y*y)
        px, py = q*q*x/psi, y/psi
        dx, dy = psi+a, psi+q*q*a
        if e < 1e-7:
            ax, ay = strength*x/dx, strength*y/dy
        else:
            ax = strength/e*np.arctan(e*x/dx)
            ay = strength/e*np.arctanh(e*y/dy)
        hxx = strength*(dx-x*px)/(dx*dx+e*e*x*x)
        hxy = -strength*x*py/(dx*dx+e*e*x*x)
        hyy = strength*(dy-y*py)/(dy*dy-e*e*y*y)
        return np.array([ax, ay, hxx, hxy, hyy])
    return nie(c["a"])-nie(c["b"])


class CombinedLens:
    def __init__(self, dark, baryons, sigma_critical):
        self.dark, self.baryons, self.sigma_critical = dark, baryons, sigma_critical

    def field(self, x, y):
        result = self.dark.field(x, y)
        for c in self.baryons:
            result += chameleon_field(x, y, c, self.sigma_critical)
        return result

    def tangential_eigenvalue(self, x, y):
        _, _, xx, xy, yy = self.field(x, y)
        return 1-(xx+yy)/2-np.hypot((xx-yy)/2, xy)

    def critical_curve(self, kpc_per_arcsec, angles=256):
        phi = np.linspace(0, 2*np.pi, angles, endpoint=False)
        scan = np.geomspace(.0001, 200, 500)
        values = self.tangential_eigenvalue(scan[:, None]*np.cos(phi), scan[:, None]*np.sin(phi))
        radii, counts = [], []
        for j, p in enumerate(phi):
            changes = np.flatnonzero((values[:-1, j] < 0) & (values[1:, j] >= 0))
            counts.append(len(changes))
            if values[-1, j] <= 0:
                raise ValueError("critical curve leaves the search domain")
            if not len(changes):
                radii.append(np.nan)
                continue
            k = changes[-1]
            radii.append(brentq(lambda r: float(self.tangential_eigenvalue(r*np.cos(p), r*np.sin(p))),
                                scan[k], scan[k+1], xtol=1e-9))
        radii = np.asarray(radii)
        all_rays = bool(np.all(np.isfinite(radii)))
        # A ray can cross multiple disjoint critical branches. Joining the
        # outermost roots across branch changes invents a closed curve and its
        # area. Preserve the envelope as a diagnostic, but do not assign it an
        # Einstein radius comparable to the published single-SIE summary.
        closed = all_rays and max(counts) == 1
        envelope = float(np.sqrt(np.mean(radii*radii))/kpc_per_arcsec) if all_rays else None
        effective = envelope if closed else None
        return {"effective_einstein_radius_arcsec": effective,
                "closed_star_shaped_outer_curve": closed,
                "maximum_outward_crossings_per_ray": max(counts),
                "outer_radial_envelope_arcsec": envelope,
                "topology_note": "single tangential branch" if closed else "multiple or missing critical branches; no single-radius comparison",
                "area_definition": "sqrt(mean(r_critical(phi)^2)); outer tangential root on every ray",
                "phi_radians": phi, "radius_kpc": radii}
