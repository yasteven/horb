// horb/orbital_basis/src/superposition.rs

/*For now, keep this intentionally boring. Don’t let it explode the model before the single-state basis works*/
use crate::{DensityField, OrbitalConfig, OrbitalError};

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
            if weight < 0.0 {
                return Err(OrbitalError::InvalidParameter {
                    name: "weight",
                    reason: "must be non-negative".into(),
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
}
