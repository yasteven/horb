// horb/orbital_basis/src/presets.rs

use crate::{OrbitalConfig, OrbitalError};

/// Milky-Way-ish 3d_z2 halo scale.
///
/// Use as a morphology starting point, not a claimed universal constant.
pub fn milky_way_dz2(dm_mass: f64) -> Result<OrbitalConfig, OrbitalError> {
    OrbitalConfig::d_z2(8.5, dm_mass)
}

/// Compact spherical test halo.
pub fn compact_1s(dm_mass: f64) -> Result<OrbitalConfig, OrbitalError> {
    OrbitalConfig::ground_state(3.0, dm_mass)
}

/// Extended spherical test halo.
pub fn extended_1s(dm_mass: f64) -> Result<OrbitalConfig, OrbitalError> {
    OrbitalConfig::ground_state(10.0, dm_mass)
}
