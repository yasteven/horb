// horb/orbital_basis/src/superposition.rs

use crate::{DensityField, OrbitalConfig, OrbitalError};

/// Signed linear density expansion: rho = sum_i weight_i M_i |psi_i|^2.
/// Weights are finite real numbers, with no sum-to-one constraint. This is
/// not the coherent quantum expression |sum_i c_i psi_i|^2: no cross terms
/// are introduced and a negative weight remains negative gravitational mass.
#[derive(Debug, Clone)]
pub struct Superposition {
    pub components: Vec<(f64, DensityField)>,
}

impl Superposition {
    pub fn new(weighted_configs: Vec<(f64, OrbitalConfig)>) -> Result<Self, OrbitalError> {
        if weighted_configs.is_empty() {
            return Err(OrbitalError::InvalidParameter {
                name: "weighted_configs",
                reason: "must contain at least one component".into(),
            });
        }

        let mut components = Vec::with_capacity(weighted_configs.len());

        for (weight, cfg) in weighted_configs {
            OrbitalConfig::new(cfg.n, cfg.l, cfg.m, cfg.a0_star, cfg.dm_mass)?;
            if !weight.is_finite() {
                return Err(OrbitalError::InvalidParameter {
                    name: "weight",
                    reason: "must be finite (positive, zero or negative)".into(),
                });
            }

            components.push((weight, DensityField::new(cfg)?));
        }

        Ok(Self { components })
    }

    pub fn rho(&self, r: f64, theta: f64, phi: f64) -> Result<f64, OrbitalError> {
        let mut total = 0.0;

        for (w, field) in &self.components {
            total += w * field.rho(r, theta, phi)?;
        }

        Ok(total)
    }

    pub fn total_mass(&self) -> f64 {
        self.components.iter().map(|(w, f)| w * f.cfg.dm_mass).sum()
    }

    /// Signed contribution to v_c^2 for the azimuthal mean of the density.
    /// Take a square root only after adding baryons and checking the total.
    pub fn disk_plane_v2(&self, r: f64, order: usize) -> Result<f64, OrbitalError> {
        self.components.iter().try_fold(0.0, |total, (w, field)| {
            Ok(total + w * crate::signed_gravity::disk_plane_v2(&field.cfg, r, order)?)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_density_and_force_are_linear_without_renormalizing() {
        let a = OrbitalConfig::ground_state(1.0, 1e11).unwrap();
        let b = OrbitalConfig::d_z2(2.0, 1e11).unwrap();
        let combined = Superposition::new(vec![(2.0, a.clone()), (-0.5, b.clone())]).unwrap();
        let fa = DensityField::new(a.clone()).unwrap();
        let fb = DensityField::new(b.clone()).unwrap();
        assert_eq!(combined.total_mass(), 1.5e11);
        for r in [1.0, 10.0, 25.0] {
            assert_eq!(
                combined.rho(r, 0.0, 0.0).unwrap(),
                2.0 * fa.rho(r, 0.0, 0.0).unwrap() - 0.5 * fb.rho(r, 0.0, 0.0).unwrap()
            );
            let ka = crate::signed_gravity::disk_plane_v2(&a, r, 96).unwrap();
            let kb = crate::signed_gravity::disk_plane_v2(&b, r, 96).unwrap();
            assert!((combined.disk_plane_v2(r, 96).unwrap() - (2.0 * ka - 0.5 * kb)).abs() < 1e-8);
        }
        let negative = Superposition::new(vec![(-1.0, a)]).unwrap();
        assert!(negative.rho(1.0, 0.0, 0.0).unwrap() < 0.0);
        assert!(negative.disk_plane_v2(1.0, 96).unwrap() < 0.0);
        assert!(Superposition::new(vec![(f64::NAN, b)]).is_err());
    }
}
