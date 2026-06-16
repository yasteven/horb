// horb/orbital_basis/src/lib.rs

//! # orbital_basis
//!
//! Hydrogenic wavefunction basis for the HORB dark-matter simulation suite.
//!
//! ## What this crate provides
//!
//! * [`OrbitalConfig`] — quantum numbers (n, l, m) plus the user-set effective
//!   Bohr radius `a0_star` (in kpc).  All physical scaling lives here so the
//!   rest of the pipeline never hard-codes units.
//!
//! * [`radial`] — associated Laguerre polynomials, radial wavefunctions R_nl,
//!   and the radial probability density P(r) = r² |R_nl|².
//!
//! * [`angular`] — real solid spherical harmonics Y_lm(θ, φ).
//!   (3d_z² = Y_2^0 is fully covered.)
//!
//! * [`density`] — full 3-D density ρ(r,θ,φ) = |ψ_nlm|², enclosed-mass
//!   tables M_enc(r), and circular-velocity profiles v_circ(R) for
//!   rotation-curve comparison.
//!
//! * [`superposition`] — weighted linear combinations Σ c_i |ψ_i|² for
//!   multi-state fits.
//!
//! * [`presets`] — named a0_star values from published galaxy–atom theories
//!   and the Fermi-bubble 3d_z² geometry match.
//!
//! ## Units convention
//!
//! | Quantity | Unit |
//! |----------|------|
//! | length   | kpc  |
//! | mass     | M_sun|
//! | velocity | km/s |
//! | G        | 4.301e-3 kpc (km/s)^2 M_sun^-1 |

/*
 * Keep M_enc(r) radial-only

Even for 3d_z2, the angular part integrates to 1 over the sphere, so the spherical enclosed mass is

M_enc(r) = M_DM ∫₀ʳ |R_nl(r')|² r'² dr'

The angular shape matters for 3D density maps / projected potentials, but the first-pass spherical rotation curve can use the radial cumulative mass.

Do not trust v_circ = sqrt(GM/R) for the final non-spherical model

That is fine for a first diagnostic, but for 3d_z2 the real circular speed should eventually come from the potential gradient in the disk plane:

v_circ²(R) = R ∂Φ(R,z=0)/∂R

So start spherical, then upgrade.
 * */
pub mod angular;
pub mod baryons;
pub mod density;
pub mod presets;
pub mod radial;
pub mod superposition;

pub use angular::RealSphericalHarmonic;
pub use baryons::BaryonicModel;
pub use density::{DensityField, EnclosedMass};
pub use radial::RadialWavefunction;
pub use superposition::Superposition;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// G in units: length = kpc, mass = M_sun, velocity = km/s
/// G = 4.3009e-6 kpc (km/s)^2 M_sun^-1
pub const G_KPC: f64 = 4.300_9e-6;

/// All quantum numbers and the effective Bohr radius for one orbital state.
///
/// `a0_star` is the *only* length scale in the problem; it sets the physical
/// size of the dark-matter halo.  Set it from a theory preset ([`presets`])
/// or supply your own value in kpc.
///
/// # Validity constraints
/// * `n >= 1`
/// * `0 <= l < n`
/// * `-l <= m <= l`
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitalConfig {
    /// Principal quantum number (n >= 1).
    pub n: u32,
    /// Azimuthal quantum number (0 <= l < n).
    pub l: u32,
    /// Magnetic quantum number (|m| <= l).
    pub m: i32,
    /// Effective Bohr radius in kpc. This is the physical scale parameter;
    /// psi(r) is evaluated at r/a0_star. Set freely to match your theory.
    pub a0_star: f64,
    /// Total dark-matter mass in M_sun. Normalises rho so M_enc(inf) = dm_mass.
    pub dm_mass: f64,
}

impl OrbitalConfig {
    /// Construct and validate quantum numbers and physical parameters.
    pub fn new(n: u32, l: u32, m: i32, a0_star: f64, dm_mass: f64) -> Result<Self, OrbitalError> {
        if n == 0 {
            return Err(OrbitalError::InvalidQuantumNumbers {
                reason: "n must be >= 1".into(),
            });
        }
        if l >= n {
            return Err(OrbitalError::InvalidQuantumNumbers {
                reason: format!("l={l} must be < n={n}"),
            });
        }
        if m.unsigned_abs() > l {
            return Err(OrbitalError::InvalidQuantumNumbers {
                reason: format!("|m|={} must be <= l={l}", m.abs()),
            });
        }
        if a0_star <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "a0_star",
                reason: "must be positive (kpc)".into(),
            });
        }
        if dm_mass <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "dm_mass",
                reason: "must be positive (M_sun)".into(),
            });
        }
        Ok(Self {
            n,
            l,
            m,
            a0_star,
            dm_mass,
        })
    }

    /// 1s ground state -- most compact, spherically symmetric DM halo.
    pub fn ground_state(a0_star: f64, dm_mass: f64) -> Result<Self, OrbitalError> {
        Self::new(1, 0, 0, a0_star, dm_mass)
    }

    /// 3d_z2 orbital: n=3, l=2, m=0.
    ///
    /// Angular part Y_2^0 ~ (3cos^2(theta) - 1), giving two polar lobes plus
    /// an equatorial torus. Morphologically matches the Fermi bubbles when
    /// a0_star is tuned to ~8-10 kpc.
    pub fn d_z2(a0_star: f64, dm_mass: f64) -> Result<Self, OrbitalError> {
        Self::new(3, 2, 0, a0_star, dm_mass)
    }

    /// Human-readable spectroscopic label, e.g. "3d_z2", "2pz", "4f+2".
    pub fn label(&self) -> String {
        let l_char = ['s', 'p', 'd', 'f', 'g', 'h']
            .get(self.l as usize)
            .copied()
            .unwrap_or('?');
        if self.m == 0 {
            match self.l {
                0 => format!("{}s", self.n),
                1 => format!("{}pz", self.n),
                2 => format!("{}d_z2", self.n),
                3 => format!("{}fz3", self.n),
                _ => format!("{}{}m0", self.n, l_char),
            }
        } else if self.m > 0 {
            format!("{}{}+{}", self.n, l_char, self.m)
        } else {
            format!("{}{}{}", self.n, l_char, self.m)
        }
    }
}

#[derive(Debug, Error)]
pub enum OrbitalError {
    #[error("invalid quantum numbers: {reason}")]
    InvalidQuantumNumbers { reason: String },

    #[error("invalid parameter '{name}': {reason}")]
    InvalidParameter { name: &'static str, reason: String },

    #[error("numerical error: {0}")]
    Numerical(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_rejects_bad_quantum_numbers() {
        assert!(OrbitalConfig::new(0, 0, 0, 1.0, 1e11).is_err());
        assert!(OrbitalConfig::new(2, 2, 0, 1.0, 1e11).is_err());
        assert!(OrbitalConfig::new(3, 2, 3, 1.0, 1e11).is_err());
    }

    #[test]
    fn config_accepts_valid_states() {
        assert!(OrbitalConfig::new(1, 0, 0, 1.0, 1e11).is_ok());
        assert!(OrbitalConfig::new(3, 2, 0, 8.5, 1e11).is_ok());
        assert!(OrbitalConfig::new(4, 3, -2, 5.0, 1e11).is_ok());
    }

    #[test]
    fn labels_are_human_readable() {
        assert_eq!(
            OrbitalConfig::ground_state(1.0, 1e11).unwrap().label(),
            "1s"
        );
        assert_eq!(OrbitalConfig::d_z2(8.5, 1e11).unwrap().label(), "3d_z2");
        assert_eq!(
            OrbitalConfig::new(2, 1, 0, 1.0, 1e11).unwrap().label(),
            "2pz"
        );
    }
}
