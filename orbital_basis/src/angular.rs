// horb/orbital_basis/src/angular.rs

use crate::OrbitalError;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct RealSphericalHarmonic {
    pub l: u32,
    pub m: i32,
}

impl RealSphericalHarmonic {
    pub fn new(l: u32, m: i32) -> Result<Self, OrbitalError> {
        if m.unsigned_abs() > l {
            return Err(OrbitalError::InvalidQuantumNumbers {
                reason: format!("|m|={} must be <= l={l}", m.abs()),
            });
        }

        Ok(Self { l, m })
    }

    /// Real spherical harmonic Y_lm(theta, phi).
    ///
    /// theta = polar angle from +z, radians.
    /// phi   = azimuthal angle in x-y plane, radians.
    pub fn value(&self, theta: f64, phi: f64) -> Result<f64, OrbitalError> {
        let c = theta.cos();
        let s = theta.sin();

        match (self.l, self.m) {
            (0, 0) => Ok(0.5 / PI.sqrt()),

            // Real p orbitals.
            (1, 0) => Ok((3.0 / (4.0 * PI)).sqrt() * c),
            (1, 1) => Ok((3.0 / (4.0 * PI)).sqrt() * s * phi.cos()),
            (1, -1) => Ok((3.0 / (4.0 * PI)).sqrt() * s * phi.sin()),

            // Real d orbitals.
            // d_z2
            (2, 0) => Ok((5.0 / (16.0 * PI)).sqrt() * (3.0 * c * c - 1.0)),

            // d_xz
            (2, 1) => Ok((15.0 / (4.0 * PI)).sqrt() * s * c * phi.cos()),

            // d_yz
            (2, -1) => Ok((15.0 / (4.0 * PI)).sqrt() * s * c * phi.sin()),

            // d_x2-y2
            (2, 2) => Ok((15.0 / (16.0 * PI)).sqrt() * s * s * (2.0 * phi).cos()),

            // d_xy
            (2, -2) => Ok((15.0 / (16.0 * PI)).sqrt() * s * s * (2.0 * phi).sin()),

            _ => Err(OrbitalError::Numerical(format!(
                "real spherical harmonic l={}, m={} not implemented yet",
                self.l, self.m
            ))),
        }
    }

    pub fn density_factor(&self, theta: f64, phi: f64) -> Result<f64, OrbitalError> {
        let y = self.value(theta, phi)?;
        Ok(y * y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn y00_value() {
        let y = RealSphericalHarmonic::new(0, 0).unwrap();
        assert_relative_eq!(y.value(0.0, 0.0).unwrap(), 0.28209479177, epsilon = 1e-10);
    }

    #[test]
    fn dz2_has_equatorial_sign_flip() {
        let y = RealSphericalHarmonic::new(2, 0).unwrap();

        let pole = y.value(0.0, 0.0).unwrap();
        let equator = y.value(std::f64::consts::FRAC_PI_2, 0.0).unwrap();

        assert!(pole > 0.0);
        assert!(equator < 0.0);
    }
}
