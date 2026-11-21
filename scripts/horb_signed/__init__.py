"""Signed linear density fitting; gravity, not speed, is the linear basis."""
from dataclasses import dataclass
import itertools
import math
import numpy as np
from scipy.optimize import minimize, linprog


@dataclass
class Basis:
    radii: np.ndarray
    matrix: np.ndarray
    metadata: list

    def at(self, radii):
        radii = np.asarray(radii)
        if not np.all(np.isfinite(radii)) or np.any(radii < self.radii[0]) or np.any(radii > self.radii[-1]):
            raise ValueError("requested radii are outside the signed basis; regenerate with a wider grid")
        return np.column_stack([np.interp(radii, self.radii, col) for col in self.matrix.T])


def load_basis(path):
    data = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8"))
    required = {"schema_version", "basis_id", "state", "n", "l", "m", "a0_star_kpc", "ref_mass_msun",
                "gravity_backend", "geometry", "normalization", "n_radial", "n_mu", "n_phi",
                "softening_kpc", "r_kpc", "v2_signed_kms2"}
    if not required.issubset(data.dtype.names or ()):
        raise ValueError("signed basis schema required (v2_signed_kms2); old speed-only libraries lost force signs and must be regenerated")
    if len(data) < 2:
        raise ValueError("empty or undersampled signed basis")
    metadata, curves, radii = [], [], None
    for bid in sorted(set(data["basis_id"])):
        rows = data[data["basis_id"] == bid]
        rows = np.sort(rows, order="r_kpc")
        meta = {}
        for key in sorted(required - {"r_kpc", "v2_signed_kms2"}):
            if not np.all(rows[key] == rows[key][0]):
                raise ValueError(f"inconsistent metadata {key} for basis {bid}")
            meta[key] = rows[key][0].item()
        for key in ("schema_version", "basis_id", "n", "l", "m", "n_radial", "n_mu", "n_phi"):
            if not np.isfinite(meta[key]) or int(meta[key]) != meta[key]:
                raise ValueError(f"{key} must be an integer")
            meta[key] = int(meta[key])
        if meta["schema_version"] != 1 or meta["geometry"] != "azimuthal_mean" or meta["normalization"] != "analytic_full_space":
            raise ValueError("unsupported signed basis convention")
        for key in ("a0_star_kpc", "ref_mass_msun"):
            if not np.isfinite(meta[key]) or meta[key] <= 0:
                raise ValueError(f"invalid {key}")
        if not (1 <= meta["n"] <= 10 and 0 <= meta["l"] < meta["n"] and meta["l"] <= 2 and abs(meta["m"]) <= meta["l"]):
            raise ValueError("unsupported quantum numbers")
        if meta["gravity_backend"] not in ("cpu_multipole", "cuda_quadrature"):
            raise ValueError("unknown gravity backend")
        if meta["basis_id"] < 0 or not (8 <= meta["n_radial"] <= 512) or not np.isfinite(meta["softening_kpc"]) or meta["softening_kpc"] < 0:
            raise ValueError("invalid basis identifier or quadrature metadata")
        if meta["gravity_backend"] == "cpu_multipole":
            if meta["n_mu"] != 8 or meta["n_phi"] != 0 or meta["softening_kpc"] != 0:
                raise ValueError("CPU multipole conventions require n_mu=8, n_phi=0, softening=0")
        elif not (8 <= meta["n_mu"] <= 512 and 8 <= meta["n_phi"] <= 512):
            raise ValueError("invalid CUDA angular quadrature metadata")
        rr, kk = rows["r_kpc"].astype(float), rows["v2_signed_kms2"].astype(float)
        if len(rr) < 2 or not np.all(np.isfinite(rr)) or np.any(rr <= 0) or np.any(np.diff(rr) <= 0) or not np.all(np.isfinite(kk)):
            raise ValueError("basis requires finite signed forces and strictly increasing positive radii")
        if radii is not None and not np.array_equal(rr, radii):
            raise ValueError("all basis curves must share the same radius grid")
        if metadata:
            for key in ("gravity_backend", "geometry", "normalization", "softening_kpc", "n_radial", "n_mu", "n_phi"):
                if meta[key] != metadata[0][key]:
                    raise ValueError(f"basis mixes {key} conventions")
        radii = rr
        metadata.append(meta)
        curves.append(kk)
    return Basis(radii, np.column_stack(curves), metadata)


def load_observations(path, lower, upper, domain):
    data = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8"))
    if not {"R_kpc", "Vobs_kms", "eV_kms"}.issubset(data.dtype.names or ()):
        raise ValueError("standard rotation-curve columns R_kpc,Vobs_kms,eV_kms are required")
    for key in ("R_kpc", "Vobs_kms", "eV_kms"):
        if not np.all(np.isfinite(data[key])):
            raise ValueError("observations must be finite")
    if np.any(data["R_kpc"] <= 0) or np.any(data["Vobs_kms"] < 0) or np.any(data["eV_kms"] <= 0):
        raise ValueError("observations require positive radii/errors and nonnegative speeds")
    if not (np.isfinite(lower) and np.isfinite(upper) and 0 < lower < upper):
        raise ValueError("invalid fit radius range")
    selected = data[(data["R_kpc"] >= lower) & (data["R_kpc"] <= upper)]
    if len(selected) < 3:
        raise ValueError("at least three observations are required in the requested fit range")
    if np.any(selected["R_kpc"] < domain[0]) or np.any(selected["R_kpc"] > domain[-1]):
        raise ValueError("fit observations exceed the basis domain; regenerate it, do not extrapolate")
    return data, np.sort(selected, order="R_kpc")


def baryons_function(path):
    rows = np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8"))
    fields = set(rows.dtype.names or ())
    if {"r_kpc", "v2_baryons_kms2"}.issubset(fields):
        # Tabulated forces support non-spherical baryons and outward gas terms.
        r, q = rows["r_kpc"].astype(float), rows["v2_baryons_kms2"].astype(float)
        if (len(r) < 2 or not np.all(np.isfinite(r)) or np.any(r <= 0)
                or np.any(np.diff(r) <= 0) or not np.all(np.isfinite(q))):
            raise ValueError("baryon force table requires finite forces and increasing positive radii")
        components = sorted(k for k in fields if k.startswith("v2_") and k.endswith("_kms2") and k != "v2_baryons_kms2")
        if components:
            parts = np.column_stack([rows[k].astype(float) for k in components])
            if not np.all(np.isfinite(parts)) or not np.allclose(parts.sum(axis=1), q, rtol=1e-10, atol=1e-8):
                raise ValueError("baryonic components do not sum to the tabulated total force")
        def tabulated(radius):
            radius = np.asarray(radius, dtype=float)
            if not np.all(np.isfinite(radius)) or np.any(radius < r[0]) or np.any(radius > r[-1]):
                raise ValueError("requested radii are outside the baryon force table; regenerate it or narrow the basis domain")
            return np.interp(radius, r, q)
        return tabulated
    if not {"kind", "mass_msun", "scale_kpc"}.issubset(rows.dtype.names or ()) or len(rows) == 0:
        raise ValueError("baryon CSV requires kind,mass_msun,scale_kpc or r_kpc,v2_baryons_kms2")
    for row in rows:
        if not np.isfinite(row["mass_msun"]) or row["mass_msun"] < 0 or not np.isfinite(row["scale_kpc"]) or row["scale_kpc"] <= 0:
            raise ValueError("baryons require nonnegative finite masses and positive finite scales")
    # The same baryonic prescription as the fixed halo comparisons.
    def evaluate(r):
        r = np.asarray(r, dtype=float)
        total = np.zeros_like(r)
        for row in rows:
            kind, mass, scale = str(row["kind"]), float(row["mass_msun"]), float(row["scale_kpc"])
            if kind in ("disk", "exponential_disk", "exponential-disk"):
                x = r / scale
                mass_enclosed = mass * (-np.expm1(-x) - x*np.exp(-x))
            elif kind in ("bulge", "hernquist", "hernquist_bulge", "hernquist-bulge"):
                mass_enclosed = mass*r*r/(r+scale)**2
            else:
                raise ValueError(f"unknown baryonic component {kind}")
            total += 4.3009e-6*mass_enclosed/r
        return total
    evaluate(np.array([1.0]))
    return evaluate


@dataclass
class Candidate:
    subset: tuple
    weights: np.ndarray
    chi2: float
    rms: float
    objective: float
    min_total_v2: float
    rank: int
    condition: float


class Search:
    """Discrete subset search with constrained continuous weight optimization.

    For nonnegative observed speeds, chi²(sqrt(b+A*w)) is convex in w on
    b+A*w>0. Weight boxes and positivity constraints are linear/convex too.
    Discrete orbital/scale selection remains a combinatorial search.
    """
    def __init__(self, A, baryons, observed, error, guard_A, guard_baryons,
                 max_weight=20.0, ridge=1e-4, criterion="bic", nonnegative=False):
        self.A, self.b, self.v, self.e = [np.asarray(x, dtype=float) for x in (A, baryons, observed, error)]
        self.C, self.d = np.asarray(guard_A, dtype=float), np.asarray(guard_baryons, dtype=float)
        if self.A.ndim != 2 or self.C.ndim != 2 or self.A.shape[1] != self.C.shape[1] or self.A.shape[0] != len(self.v):
            raise ValueError("inconsistent fit matrices")
        if self.b.shape != self.v.shape or self.e.shape != self.v.shape or self.d.shape != (self.C.shape[0],):
            raise ValueError("inconsistent observation/baryon arrays")
        if not all(np.all(np.isfinite(x)) for x in (self.A,self.b,self.v,self.e,self.C,self.d)) or np.any(self.e<=0) or np.any(self.v<0):
            raise ValueError("invalid fit data")
        if not np.isfinite(max_weight) or max_weight<=0 or not np.isfinite(ridge) or ridge<0 or criterion not in ("bic","chi2"):
            raise ValueError("invalid bounds, ridge or selection criterion")
        self.bound, self.ridge, self.criterion, self.nonnegative = max_weight, ridge, criterion, nonnegative
        self.cache = {}
        self.failed = 0

    def fit(self, subset):
        subset = tuple(sorted(subset))
        if subset in self.cache:
            return self.cache[subset]
        k, n = len(subset), len(self.v)
        A, C = self.A[:,subset], self.C[:,subset]
        if k == 0:
            if np.any(self.b<0) or np.any(self.d<0):
                self.cache[subset] = None
                return None
            w, q = np.empty(0), self.b
        else:
            bounds = [(0 if self.nonnegative else -self.bound, self.bound)]*k
            # Squared-speed solution is an initializer only. Negative dark
            # targets b>vobs² are preserved, and final fitting is in velocity.
            sigma_q = np.maximum(2*self.v*self.e, self.e**2)
            w0 = np.linalg.lstsq(A/sigma_q[:,None], (self.v**2-self.b)/sigma_q, rcond=1e-12)[0]
            w0 = np.clip(w0, bounds[0][0], self.bound)
            floor = 1e-8
            if np.min(self.d+C@w0)<floor:
                if np.min(self.d)>=floor:
                    w0 = np.zeros(k)
                else:
                    feasible = linprog(np.zeros(k), A_ub=-C/1e4, b_ub=(self.d-floor)/1e4, bounds=bounds, method="highs")
                    if not feasible.success:
                        self.failed += 1
                        self.cache[subset] = None
                        return None
                    w0 = feasible.x
            def objective(w):
                q = self.b+A@w
                root = np.sqrt(np.maximum(q, floor))
                dq = (1-self.v/root)/self.e**2
                # Linear continuation below the feasible domain is only for
                # optimizer trial steps. Such points can never be accepted.
                val = np.sum(((root-self.v)/self.e)**2 + dq*np.minimum(q-floor,0)) + self.ridge*(w@w)
                grad = A.T@dq + 2*self.ridge*w
                return val/n, grad/n
            constraint = {"type":"ineq", "fun":lambda w:(self.d+C@w-floor)/1e4,
                          "jac":lambda w:C/1e4}
            result = minimize(objective, w0, jac=True, method="SLSQP", bounds=bounds,
                              constraints=[constraint], options={"ftol":1e-10,"maxiter":500})
            w = result.x
            q = self.b+A@w
            if (not result.success or not np.all(np.isfinite(w)) or np.any(q<0)
                    or np.any(self.d+C@w<0) or np.any(w<bounds[0][0]-1e-8) or np.any(w>self.bound+1e-8)):
                self.failed += 1
                self.cache[subset] = None
                return None
        residual = np.sqrt(q)-self.v
        chi2 = float(np.sum((residual/self.e)**2))
        score = chi2+self.ridge*float(w@w) + (k*math.log(n) if self.criterion=="bic" else 0)
        rank = int(np.linalg.matrix_rank(A/self.e[:,None])) if k else 0
        condition = float(np.linalg.cond(A/self.e[:,None])) if k else 1.0
        candidate = Candidate(subset,w,chi2,float(np.sqrt(np.mean(residual**2))),score,
                              float(np.min(self.d+C@w)),rank,condition)
        self.cache[subset] = candidate
        return candidate

    def run(self, max_components=4, method="beam", beam_width=6, budget=2000, seed=42):
        count = self.A.shape[1]
        if max_components<0 or max_components>count or beam_width<1 or budget<1 or method not in ("beam","anneal","exhaustive"):
            raise ValueError("invalid subset search controls")
        if budget<count+1 and max_components>0:
            raise ValueError("budget must accommodate the null model and every single orbital")
        if method=="exhaustive" and sum(math.comb(count,k) for k in range(max_components+1))>budget:
            raise ValueError("exhaustive combinations exceed budget; increase budget or use beam/anneal")
        self.fit(())
        if max_components:
            for i in range(count): self.fit((i,))
        if method=="exhaustive":
            for k in range(2,max_components+1):
                for subset in itertools.combinations(range(count),k): self.fit(subset)
        elif method=="beam":
            frontier = [(i,) for i in range(count)] if max_components else []
            for k in range(2,max_components+1):
                frontier = sorted(frontier,key=lambda s:(self.cache[s].objective if self.cache[s] else math.inf,s))[:beam_width]
                proposals = sorted({tuple(sorted((*parent,i))) for parent in frontier for i in range(count) if i not in parent})
                frontier = []
                for subset in proposals:
                    if len(self.cache)>=budget: break
                    self.fit(subset)
                    frontier.append(subset)
        elif max_components:
            rng = np.random.default_rng(seed)
            valid = [r for r in self.cache.values() if r is not None]
            current = min(valid,key=lambda r:(r.objective,r.subset)) if valid else None
            # Add/drop/swap proposals choose orbitals and their count; each
            # proposal gets a full constrained weight refit. Seed is recorded.
            for step in range(budget*10):
                if len(self.cache)>=budget: break
                if current is None or step%10==0:
                    size = int(rng.integers(1,max_components+1))
                    proposal = tuple(sorted(rng.choice(count,size,replace=False).tolist()))
                else:
                    proposed = set(current.subset)
                    move = int(rng.integers(3))
                    if proposed and move in (0,2): proposed.remove(int(rng.choice(sorted(proposed))))
                    remaining = sorted(set(range(count))-proposed)
                    if remaining and len(proposed)<max_components and move in (1,2): proposed.add(int(rng.choice(remaining)))
                    proposal = tuple(sorted(proposed))
                candidate = self.fit(proposal)
                temperature = max(0.01,10*(1-len(self.cache)/budget)**2)
                if candidate and (current is None or candidate.objective<current.objective or rng.random()<math.exp(min(0,(current.objective-candidate.objective)/temperature))):
                    current = candidate
        valid = sorted((r for r in self.cache.values() if r is not None),key=lambda r:(r.objective,r.subset))
        if not valid: raise ValueError("no converged feasible fit; inspect basis, baryons and coefficient bounds")
        return valid
