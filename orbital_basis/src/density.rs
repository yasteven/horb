// horb/orbital_basis/src/density.rs

use crate::{
    angular::RealSphericalHarmonic, radial::RadialWavefunction, OrbitalConfig, OrbitalError, G_KPC,
};

#[derive(Debug, Clone)]
pub struct DensityField {
    pub cfg: OrbitalConfig,
    radial: RadialWavefunction,
    angular: RealSphericalHarmonic,
}

impl DensityField {
    pub fn new(cfg: OrbitalConfig) -> Result<Self, OrbitalError> {
        let radial = RadialWavefunction::new(cfg.clone());
        let angular = RealSphericalHarmonic::new(cfg.l, cfg.m)?;

        Ok(Self {
            cfg,
            radial,
            angular,
        })
    }

    /// Full 3D density:
    ///
    /// rho(r,theta,phi) = M_DM |R_nl(r)|² |Y_lm(theta,phi)|²
    ///
    /// Units: M_sun / kpc³
    pub fn rho(&self, r: f64, theta: f64, phi: f64) -> Result<f64, OrbitalError> {
        let radial = self.radial.radial_density_factor(r)?;
        let angular = self.angular.density_factor(theta, phi)?;

        Ok(self.cfg.dm_mass * radial * angular)
    }

    /// Spherical enclosed mass by radial integration.
    ///
    /// Since ∫ |Y_lm|² dΩ = 1, only the radial probability matters.
    pub fn enclosed_mass(&self, r_max: f64, steps: usize) -> Result<f64, OrbitalError> {
        if r_max < 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "r_max",
                reason: "must be non-negative".into(),
            });
        }

        if steps < 2 {
            return Err(OrbitalError::InvalidParameter {
                name: "steps",
                reason: "must be >= 2".into(),
            });
        }

        let h = r_max / steps as f64;
        let mut sum = 0.0;

        for i in 0..=steps {
            let r = i as f64 * h;
            let weight = if i == 0 || i == steps {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };

            sum += weight * self.radial.radial_probability_density(r)?;
        }

        let integral = sum * h / 3.0;

        Ok(self.cfg.dm_mass * integral)
    }

    /// First-pass spherical circular velocity:
    ///
    /// v_circ(r) = sqrt(G M_enc(r) / r)
    ///
    /// Good for debugging and spherical comparison.
    /// Not the final non-spherical disk-plane velocity.
    pub fn circular_velocity_spherical(&self, r: f64, steps: usize) -> Result<f64, OrbitalError> {
        if r <= 0.0 {
            return Ok(0.0);
        }

        let m = self.enclosed_mass(r, steps)?;
        Ok((G_KPC * m / r).sqrt())
    }
}

#[derive(Debug, Clone)]
pub struct EnclosedMass {
    pub radii_kpc: Vec<f64>,
    pub mass_msun: Vec<f64>,
}

impl EnclosedMass {
    pub fn build(
        field: &DensityField,
        r_max: f64,
        bins: usize,
        integrate_steps: usize,
    ) -> Result<Self, OrbitalError> {
        if bins < 2 {
            return Err(OrbitalError::InvalidParameter {
                name: "bins",
                reason: "must be >= 2".into(),
            });
        }

        let mut radii_kpc = Vec::with_capacity(bins);
        let mut mass_msun = Vec::with_capacity(bins);

        for i in 0..bins {
            let r = r_max * i as f64 / (bins - 1) as f64;
            radii_kpc.push(r);
            mass_msun.push(field.enclosed_mass(r, integrate_steps)?);
        }

        Ok(Self {
            radii_kpc,
            mass_msun,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn density_is_positive() {
        let cfg = OrbitalConfig::d_z2(8.5, 1e11).unwrap();
        let field = DensityField::new(cfg).unwrap();

        let rho = field.rho(8.0, 0.1, 0.0).unwrap();
        assert!(rho >= 0.0);
        assert!(rho.is_finite());
    }

    #[test]
    fn enclosed_mass_approaches_total_mass() {
        let cfg = OrbitalConfig::ground_state(1.0, 1e11).unwrap();
        let field = DensityField::new(cfg).unwrap();

        let m = field.enclosed_mass(30.0, 10_000).unwrap();

        assert!(m > 0.999 * 1e11);
        assert!(m < 1.001 * 1e11);
    }

    #[test]
    fn circular_velocity_is_finite() {
        let cfg = OrbitalConfig::ground_state(3.0, 1e11).unwrap();
        let field = DensityField::new(cfg).unwrap();

        let v = field.circular_velocity_spherical(8.0, 2000).unwrap();

        assert!(v.is_finite());
        assert!(v > 0.0);
    }
}
