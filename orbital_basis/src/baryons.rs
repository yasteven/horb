// horb/orbital_basis/src/baryons.rs

use crate::{OrbitalError, G_KPC};

/// Simple baryonic placeholder model.
///
/// This is intentionally not final galaxy dynamics.
/// It gives the pipeline a clean v_total² = v_dm² + v_disk² + v_bulge² path.
///
/// Disk: sphericalized exponential enclosed-mass approximation
///   M_disk(<R) = M_disk [1 - exp(-R/Rd)(1 + R/Rd)]
///
/// Bulge: Hernquist enclosed mass
///   M_bulge(<R) = M_bulge R² / (R + a)²
#[derive(Debug, Clone, Copy)]
pub struct BaryonicModel {
    pub disk_mass_msun: f64,
    pub disk_scale_kpc: f64,
    pub bulge_mass_msun: f64,
    pub bulge_scale_kpc: f64,
}

impl BaryonicModel {
    pub fn new(
        disk_mass_msun: f64,
        disk_scale_kpc: f64,
        bulge_mass_msun: f64,
        bulge_scale_kpc: f64,
    ) -> Result<Self, OrbitalError> {
        if disk_mass_msun < 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "disk_mass_msun",
                reason: "must be non-negative".into(),
            });
        }

        if bulge_mass_msun < 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "bulge_mass_msun",
                reason: "must be non-negative".into(),
            });
        }

        if disk_mass_msun > 0.0 && disk_scale_kpc <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "disk_scale_kpc",
                reason: "must be positive when disk mass is nonzero".into(),
            });
        }

        if bulge_mass_msun > 0.0 && bulge_scale_kpc <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "bulge_scale_kpc",
                reason: "must be positive when bulge mass is nonzero".into(),
            });
        }

        Ok(Self {
            disk_mass_msun,
            disk_scale_kpc,
            bulge_mass_msun,
            bulge_scale_kpc,
        })
    }

    pub fn disk_enclosed_mass(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 || self.disk_mass_msun == 0.0 {
            return 0.0;
        }

        let x = r_kpc / self.disk_scale_kpc;
        self.disk_mass_msun * (1.0 - (-x).exp() * (1.0 + x))
    }

    pub fn bulge_enclosed_mass(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 || self.bulge_mass_msun == 0.0 {
            return 0.0;
        }

        self.bulge_mass_msun * r_kpc.powi(2) / (r_kpc + self.bulge_scale_kpc).powi(2)
    }

    pub fn disk_velocity(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        let m = self.disk_enclosed_mass(r_kpc);
        (G_KPC * m / r_kpc).sqrt()
    }

    pub fn bulge_velocity(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        let m = self.bulge_enclosed_mass(r_kpc);
        (G_KPC * m / r_kpc).sqrt()
    }

    pub fn baryon_velocity(&self, r_kpc: f64) -> f64 {
        let vd = self.disk_velocity(r_kpc);
        let vb = self.bulge_velocity(r_kpc);

        (vd * vd + vb * vb).sqrt()
    }

    pub fn total_velocity(
        dm_velocity_kms: f64,
        disk_velocity_kms: f64,
        bulge_velocity_kms: f64,
    ) -> f64 {
        (dm_velocity_kms.powi(2) + disk_velocity_kms.powi(2) + bulge_velocity_kms.powi(2)).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn baryonic_model_accepts_reasonable_values() {
        let b = BaryonicModel::new(6e10, 3.0, 1e10, 0.7).unwrap();

        assert!(b.disk_velocity(8.0).is_finite());
        assert!(b.bulge_velocity(8.0).is_finite());
        assert!(b.baryon_velocity(8.0) > 0.0);
    }

    #[test]
    fn zero_baryons_give_zero_velocity() {
        let b = BaryonicModel::new(0.0, 0.0, 0.0, 0.0).unwrap();

        assert_eq!(b.disk_velocity(8.0), 0.0);
        assert_eq!(b.bulge_velocity(8.0), 0.0);
        assert_eq!(b.baryon_velocity(8.0), 0.0);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaryonComponentKind {
    ExponentialDisk,
    HernquistBulge,
}

impl BaryonComponentKind {
    pub fn parse(value: &str) -> Option<Self> {
        match value.trim().to_lowercase().as_str() {
            "disk" | "exponential_disk" | "exponential-disk" => Some(Self::ExponentialDisk),
            "bulge" | "hernquist" | "hernquist_bulge" | "hernquist-bulge" => {
                Some(Self::HernquistBulge)
            }
            _ => None,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            Self::ExponentialDisk => "exponential_disk",
            Self::HernquistBulge => "hernquist_bulge",
        }
    }
}

#[derive(Debug, Clone)]
pub struct BaryonComponent {
    pub label: String,
    pub kind: BaryonComponentKind,
    pub mass_msun: f64,
    pub scale_kpc: f64,
}

impl BaryonComponent {
    pub fn new(
        label: impl Into<String>,
        kind: BaryonComponentKind,
        mass_msun: f64,
        scale_kpc: f64,
    ) -> anyhow::Result<Self> {
        if mass_msun < 0.0 {
            anyhow::bail!("baryon component mass must be non-negative");
        }

        if scale_kpc <= 0.0 {
            anyhow::bail!("baryon component scale must be positive");
        }

        Ok(Self {
            label: label.into(),
            kind,
            mass_msun,
            scale_kpc,
        })
    }

    pub fn enclosed_mass(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        match self.kind {
            BaryonComponentKind::ExponentialDisk => {
                let x = r_kpc / self.scale_kpc;
                self.mass_msun * (1.0 - (-x).exp() * (1.0 + x))
            }
            BaryonComponentKind::HernquistBulge => {
                self.mass_msun * r_kpc.powi(2) / (r_kpc + self.scale_kpc).powi(2)
            }
        }
    }

    pub fn velocity(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        (crate::G_KPC * self.enclosed_mass(r_kpc) / r_kpc).sqrt()
    }
}

#[derive(Debug, Clone)]
pub struct MultiBaryonicModel {
    pub components: Vec<BaryonComponent>,
}

impl MultiBaryonicModel {
    pub fn new(components: Vec<BaryonComponent>) -> anyhow::Result<Self> {
        if components.is_empty() {
            anyhow::bail!("multi-baryonic model requires at least one component");
        }

        Ok(Self { components })
    }

    pub fn enclosed_mass(&self, r_kpc: f64) -> f64 {
        self.components
            .iter()
            .map(|component| component.enclosed_mass(r_kpc))
            .sum()
    }

    pub fn baryon_velocity(&self, r_kpc: f64) -> f64 {
        let v2: f64 = self
            .components
            .iter()
            .map(|component| component.velocity(r_kpc).powi(2))
            .sum();

        v2.sqrt()
    }

    pub fn total_velocity(&self, r_kpc: f64, v_dm_kms: f64) -> f64 {
        let v_baryon = self.baryon_velocity(r_kpc);
        (v_dm_kms * v_dm_kms + v_baryon * v_baryon).sqrt()
    }

    pub fn component_velocities(&self, r_kpc: f64) -> Vec<(String, f64)> {
        self.components
            .iter()
            .map(|component| (component.label.clone(), component.velocity(r_kpc)))
            .collect()
    }
}
