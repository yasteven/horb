"""McCulloch (2017), DOI 10.1007/s10509-017-3128-6, Equations 8 and 9.

The published acceleration prescription gives the total speed from baryons.
It is not an additional dark halo and supplies no lensing prescription here.
"""
import numpy as np

PAPER_URL = "https://doi.org/10.1007/s10509-017-3128-6"
C_M_S = 299792458.0
THETA_M = 8.8e26
KPC_M = 3.085677581491367e19
A_QI_M_S2 = 2 * C_M_S**2 / THETA_M


def acceleration_scale(redshift=0.0):
    """Use the paper's Eq. 9, Theta(z) = Theta(0)/(1+z)."""
    z = float(redshift)
    if not np.isfinite(z) or z < 0:
        raise ValueError("QI requires a finite nonnegative cosmological redshift")
    return A_QI_M_S2 * (1 + z)


def velocity(r_kpc, baryons_v2_kms2, redshift=0.0):
    """Total circular speed [km/s]; radius [kpc], inward baryonic v² [(km/s)²]."""
    r, b = np.broadcast_arrays(np.asarray(r_kpc, dtype=float),
                               np.asarray(baryons_v2_kms2, dtype=float))
    if not np.isfinite(r).all() or np.any(r <= 0):
        raise ValueError("QI requires finite positive radii")
    if not np.isfinite(b).all() or np.any(b < 0):
        raise ValueError("QI Equation 8 requires finite nonnegative total baryonic v^2")
    # v^4 = v_bar^4 + a_QI*r*v_bar^2. Factored roots also preserve the
    # continuous zero-force limit without a division by zero.
    a = acceleration_scale(redshift)
    return np.sqrt(np.sqrt(b) * np.sqrt(b + a * KPC_M / 1e6 * r))


def metadata(redshift=0.0):
    a = acceleration_scale(redshift)
    return {"source": PAPER_URL, "equation": 8, "redshift_scaling_equation": 9,
            "redshift": float(redshift), "c_m_s": C_M_S,
            "theta_now_m": THETA_M, "theta_m": THETA_M / (1 + float(redshift)),
            "a_qi_m_s2": a, "parameter_count": 0, "fitted": False,
            "baryons": "Same fixed baryonic force; QI is the total speed, with no additional halo.",
            "domain": "Positive radius and nonnegative total inward baryonic v^2."}
