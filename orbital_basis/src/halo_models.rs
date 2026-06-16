// horb/orbital_basis/src/halo_models.rs

use crate::{OrbitalError, G_KPC};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClassicalHaloKind {
    PseudoIsothermal,
    Nfw,
    Burkert,
}

impl ClassicalHaloKind {
    pub fn parse(s: &str) -> Result<Self, OrbitalError> {
        match s {
            "piso" | "pseudo_isothermal" | "pseudo-isothermal" => Ok(Self::PseudoIsothermal),
            "nfw" | "NFW" => Ok(Self::Nfw),
            "burkert" => Ok(Self::Burkert),
            _ => Err(OrbitalError::InvalidParameter {
                name: "classical_halo_kind",
                reason: format!("unknown halo kind '{s}'; use piso, nfw, or burkert"),
            }),
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            Self::PseudoIsothermal => "piso",
            Self::Nfw => "nfw",
            Self::Burkert => "burkert",
        }
    }
}

/// Classical spherical halo profile normalized by enclosed mass.
///
/// Parameters:
///
/// scale_kpc:
///   core radius for pseudo-isothermal,
///   scale radius for NFW,
///   core radius for Burkert.
///
/// density_norm:
///   rho0 for pseudo-isothermal and Burkert,
///   rho_s for NFW.
///
/// The constructor `from_mass_at_radius` solves density_norm so that:
///
///   M(<r_ref_kpc) = m_ref_msun
#[derive(Debug, Clone, Copy)]
pub struct ClassicalHalo {
    pub kind: ClassicalHaloKind,
    pub scale_kpc: f64,
    pub density_norm_msun_per_kpc3: f64,
}

impl ClassicalHalo {
    pub fn from_mass_at_radius(
        kind: ClassicalHaloKind,
        scale_kpc: f64,
        m_ref_msun: f64,
        r_ref_kpc: f64,
    ) -> Result<Self, OrbitalError> {
        if scale_kpc <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "scale_kpc",
                reason: "must be positive".into(),
            });
        }

        if m_ref_msun <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "m_ref_msun",
                reason: "must be positive".into(),
            });
        }

        if r_ref_kpc <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "r_ref_kpc",
                reason: "must be positive".into(),
            });
        }

        let unit_mass = Self::unit_density_mass(kind, scale_kpc, r_ref_kpc);

        if unit_mass <= 0.0 || !unit_mass.is_finite() {
            return Err(OrbitalError::Numerical(
                "failed to compute positive finite unit-density halo mass".into(),
            ));
        }

        let density_norm_msun_per_kpc3 = m_ref_msun / unit_mass;

        Ok(Self {
            kind,
            scale_kpc,
            density_norm_msun_per_kpc3,
        })
    }

    pub fn density(&self, r_kpc: f64) -> f64 {
        if r_kpc < 0.0 {
            return f64::NAN;
        }

        let r = r_kpc.max(1.0e-12);
        let x = r / self.scale_kpc;
        let rho = self.density_norm_msun_per_kpc3;

        match self.kind {
            ClassicalHaloKind::PseudoIsothermal => rho / (1.0 + x * x),
            ClassicalHaloKind::Nfw => rho / (x * (1.0 + x).powi(2)),
            ClassicalHaloKind::Burkert => rho / ((1.0 + x) * (1.0 + x * x)),
        }
    }

    pub fn enclosed_mass(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        self.density_norm_msun_per_kpc3 * Self::unit_density_mass(self.kind, self.scale_kpc, r_kpc)
    }

    pub fn circular_velocity(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        let m = self.enclosed_mass(r_kpc);
        (G_KPC * m / r_kpc).sqrt()
    }

    fn unit_density_mass(kind: ClassicalHaloKind, scale_kpc: f64, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        let a = scale_kpc;
        let x = r_kpc / a;

        match kind {
            ClassicalHaloKind::PseudoIsothermal => {
                // rho = rho0 / (1 + x^2)
                // M(r) = 4 pi rho0 a^3 [x - atan(x)]
                4.0 * PI * a.powi(3) * (x - x.atan())
            }

            ClassicalHaloKind::Nfw => {
                // rho = rho_s / [x(1+x)^2]
                // M(r) = 4 pi rho_s a^3 [ln(1+x) - x/(1+x)]
                4.0 * PI * a.powi(3) * ((1.0 + x).ln() - x / (1.0 + x))
            }

            ClassicalHaloKind::Burkert => {
                // rho = rho0 / [(1+x)(1+x^2)]
                // M(r) = pi rho0 a^3 [2ln(1+x) + ln(1+x^2) - 2atan(x)]
                PI * a.powi(3) * (2.0 * (1.0 + x).ln() + (1.0 + x * x).ln() - 2.0 * x.atan())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classical_halos_normalize_to_reference_mass() {
        for kind in [
            ClassicalHaloKind::PseudoIsothermal,
            ClassicalHaloKind::Nfw,
            ClassicalHaloKind::Burkert,
        ] {
            let h = ClassicalHalo::from_mass_at_radius(kind, 10.0, 1e11, 80.0).unwrap();
            let m = h.enclosed_mass(80.0);

            let rel = (m - 1e11).abs() / 1e11;
            assert!(rel < 1e-12, "kind={:?}, rel={}", kind, rel);
        }
    }

    #[test]
    fn classical_halo_velocities_are_finite() {
        let h = ClassicalHalo::from_mass_at_radius(ClassicalHaloKind::Burkert, 10.0, 1e11, 80.0)
            .unwrap();

        let v = h.circular_velocity(8.0);

        assert!(v.is_finite());
        assert!(v > 0.0);
    }
}
