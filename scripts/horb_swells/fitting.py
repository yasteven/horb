"""Fit rotation measurements only; this module never reads lensing outcomes."""
import itertools
import numpy as np
from scipy.optimize import differential_evolution, least_squares, minimize
from horb_highlights import metrics
from .models import baseline_v2


def fit_classical(model, grid, baryons, projection, observed, error):
    bounds = [(8., 13.), (np.log10(.05), np.log10(50 if model == "single" else 500))]
    if model == "einasto":
        bounds.append((np.log10(.03), np.log10(2.)))
    def residual(p):
        speed = np.sqrt(np.maximum(baryons+baseline_v2(model, grid, p), 0))
        return (projection@speed-observed)/error
    global_fit = differential_evolution(lambda p: residual(p)@residual(p), bounds,
                                        seed=42, popsize=12, maxiter=400, tol=1e-8, polish=False)
    local = least_squares(residual, global_fit.x, bounds=np.array(bounds).T,
                          ftol=1e-11, xtol=1e-11, gtol=1e-11, max_nfev=3000)
    if not local.success:
        raise RuntimeError(f"{model} optimization failed: {local.message}")
    prediction = observed+local.fun*error
    return {**metrics(prediction, observed, error), "parameters": (10**local.x).tolist(),
            "log10_parameters": local.x.tolist(), "log10_bounds": bounds,
            "mass_convention": "total_mass" if model == "single" else "mass_within_13_kpc",
            "reference_radius_kpc": 13., "parameter_count": len(bounds),
            "bound_reached": any(min(x-a, b-x) < 1e-4 for x, (a, b) in zip(local.x, bounds)),
            "global_converged": bool(global_fit.success), "local_converged": bool(local.success)}


class ProjectedSearch:
    """Subset beam search with bounded signed weights and inward total force.

    PSF/slit integration is nonlinear in force coefficients, so the old
    direct-circular-speed convexity argument does not apply. Use two feasible
    starts for each subset and report the best converged searched candidate.
    No lensing radius, ellipticity, or image is used for selection.
    """
    def __init__(self, matrix, baryons, projection, observed, error,
                 max_weight=20., nonnegative=False, validation_matrix=None, validation_baryons=None):
        self.A, self.b, self.P = matrix, baryons, projection
        self.v, self.e = observed, error
        self.bound, self.nonnegative = max_weight, nonnegative
        self.cache, self.failed = {}, 0
        # Origin has exactly zero force for all models. All other radii guard
        # against negative total v², including regions outside the fitted data.
        self.guard = np.arange(1, len(baryons))
        self.validation_A = matrix if validation_matrix is None else validation_matrix
        self.validation_b = baryons if validation_baryons is None else validation_baryons

    def fit(self, subset):
        subset = tuple(sorted(subset))
        if subset in self.cache:
            return self.cache[subset]
        if not subset:
            result = {**metrics(self.P@np.sqrt(self.b), self.v, self.e),
                      "subset": [], "weights": [], "bound_reached": False}
            self.cache[subset] = result
            return result
        A = np.ascontiguousarray(self.A[:, subset])
        guard = self.guard
        C, d = A[guard], self.b[guard]
        validation_A = self.validation_A[:, subset]
        validation_b = self.validation_b
        # Normalize constraints row by row to handle the tiny central force.
        normal = np.maximum(np.abs(d)+np.max(np.abs(C), axis=1), 1e-12)
        Cn, dn = C/normal[:, None], d/normal
        lower = 0 if self.nonnegative else -self.bound
        def objective(w):
            q = self.b+A@w
            root = np.sqrt(np.maximum(q, 1e-20))
            residual = (self.P@root-self.v)/self.e
            jac = self.P@(A/(2*root[:, None]))
            return float(residual@residual), 2*jac.T@(residual/self.e)
        # Linearize the complete observing model at the baryon curve, then
        # shrink this initial step until it obeys every force constraint.
        root = np.sqrt(np.maximum(self.b, 1e-20))
        jac = self.P@(A/(2*root[:, None]))
        w0 = np.linalg.lstsq(jac/self.e[:, None], (self.v-self.P@root)/self.e, rcond=1e-10)[0]
        w0 = np.clip(w0, lower, self.bound)
        for _ in range(80):
            if np.min(dn+Cn@w0) >= 0:
                break
            w0 *= .5
        candidates = []
        for start in (w0, np.zeros(len(subset))):
            # Refine constraints where cancellation causes a force minimum
            # between the coarser observing-grid nodes. Only kinematics enter.
            for _ in range(8):
                fit = minimize(objective, start, jac=True, method="SLSQP",
                               bounds=[(lower, self.bound)]*len(subset),
                               constraints=[{"type": "ineq", "fun": lambda w: dn+Cn@w,
                                             "jac": lambda w: Cn}],
                               options={"maxiter": 250, "ftol": 1e-9})
                fine_force = validation_b+validation_A@fit.x
                bad = np.flatnonzero(fine_force < -1e-7)
                if not fit.success or not len(bad):
                    break
                extra_A, extra_b = validation_A[bad], validation_b[bad]
                extra_norm = np.maximum(np.abs(extra_b)+np.max(np.abs(extra_A), axis=1), 1e-12)
                Cn = np.vstack([Cn, extra_A/extra_norm[:, None]])
                # A tiny positive buffer avoids repeatedly accepting numerical
                # negative zero in SLSQP's absolute feasibility tolerance.
                dn = np.r_[dn, (extra_b-1e-5)/extra_norm]
                start = fit.x
            if (fit.success and np.all(np.isfinite(fit.x))
                    and np.min(self.b+A@fit.x) >= -1e-7
                    and np.min(fine_force) >= -1e-7
                    and np.min(dn+Cn@fit.x) >= -1e-10):
                prediction = self.P@np.sqrt(np.maximum(self.b+A@fit.x, 0))
                candidates.append({**metrics(prediction, self.v, self.e), "subset": list(subset),
                                   "weights": fit.x.tolist(),
                                   "bound_reached": bool(np.any(fit.x-lower < 1e-4)
                                                         or np.any(self.bound-fit.x < 1e-4))})
        result = min(candidates, key=lambda c: c["chi2"]) if candidates else None
        self.cache[subset] = result
        if result is None:
            self.failed += 1
        return result

    def run(self, max_components=4, beam_width=8, budget=1200, exhaustive=False):
        size = self.A.shape[1]
        if (not 0 <= max_components <= min(4, size) or beam_width < 1
                or budget < (1 if max_components == 0 else size+1)):
            raise ValueError("invalid subset search controls")
        if exhaustive:
            import math
            if sum(math.comb(size, k) for k in range(max_components+1)) > budget:
                raise ValueError("budget too small for exhaustive subset search")
        self.fit(())
        if max_components == 0:
            return [self.cache[()]]
        frontier = []
        for i in range(size):
            if self.fit((i,)) is not None:
                frontier.append((i,))
        for count in range(2, max_components+1):
            if exhaustive:
                proposals = itertools.combinations(range(size), count)
            else:
                frontier.sort(key=lambda s: (self.cache[s]["chi2"], s))
                proposals = sorted({tuple(sorted((*parent, i))) for parent in frontier[:beam_width]
                                    for i in range(size) if i not in parent})
            frontier = []
            for subset in proposals:
                if len(self.cache) >= budget:
                    break
                if self.fit(subset) is not None:
                    frontier.append(subset)
        return sorted((c for c in self.cache.values() if c is not None),
                      key=lambda c: (c["chi2"], c["subset"]))
