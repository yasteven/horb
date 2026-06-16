// horb/orbital_basis/src/radial.rs

use crate::{OrbitalConfig, OrbitalError};

#[derive(Debug, Clone)]
pub struct RadialWavefunction {
    pub cfg: OrbitalConfig,
}

impl RadialWavefunction {
    pub fn new(cfg: OrbitalConfig) -> Self {
        Self { cfg }
    }

    /// Hydrogenic radial wavefunction R_nl(r), normalized so:
    ///
    /// ∫₀∞ |R_nl(r)|² r² dr = 1
    ///
    /// r and a0_star are in kpc, so R has units kpc^(-3/2).
    pub fn value(&self, r: f64) -> Result<f64, OrbitalError> {
        if r < 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "r",
                reason: "radius must be non-negative".into(),
            });
        }

        let n = self.cfg.n as i32;
        let l = self.cfg.l as i32;
        let a = self.cfg.a0_star;

        let k = n - l - 1;
        let alpha = 2 * l + 1;

        let rho = 2.0 * r / ((n as f64) * a);

        let norm = ((2.0 / ((n as f64) * a)).powi(3) * factorial(k)? as f64
            / (2.0 * n as f64 * factorial(n + l)? as f64))
            .sqrt();

        let lag = associated_laguerre(k, alpha, rho)?;

        Ok(norm * (-rho / 2.0).exp() * rho.powi(l) * lag)
    }

    /// Radial probability density:
    ///
    /// P(r) = r² |R_nl(r)|²
    ///
    /// This integrates to 1 over dr.
    pub fn radial_probability_density(&self, r: f64) -> Result<f64, OrbitalError> {
        let rnl = self.value(r)?;
        Ok(r * r * rnl * rnl)
    }

    /// Raw spatial radial density factor |R_nl(r)|².
    pub fn radial_density_factor(&self, r: f64) -> Result<f64, OrbitalError> {
        let rnl = self.value(r)?;
        Ok(rnl * rnl)
    }
}

/// Associated Laguerre polynomial L_k^alpha(x).
///
/// Uses recurrence:
///
/// L_0^a = 1
/// L_1^a = 1 + a - x
///
/// k L_k^a = (2k - 1 + a - x)L_{k-1}^a - (k - 1 + a)L_{k-2}^a
pub fn associated_laguerre(k: i32, alpha: i32, x: f64) -> Result<f64, OrbitalError> {
    if k < 0 {
        return Err(OrbitalError::Numerical(
            "associated Laguerre requires k >= 0".into(),
        ));
    }

    if k == 0 {
        return Ok(1.0);
    }

    if k == 1 {
        return Ok(1.0 + alpha as f64 - x);
    }

    let mut lm2 = 1.0;
    let mut lm1 = 1.0 + alpha as f64 - x;
    let mut lcur = lm1;

    for j in 2..=k {
        let jf = j as f64;
        lcur = ((2.0 * jf - 1.0 + alpha as f64 - x) * lm1 - (jf - 1.0 + alpha as f64) * lm2) / jf;

        lm2 = lm1;
        lm1 = lcur;
    }

    Ok(lcur)
}

fn factorial(n: i32) -> Result<u64, OrbitalError> {
    if n < 0 {
        return Err(OrbitalError::Numerical(
            "factorial called with negative input".into(),
        ));
    }

    let mut acc: u64 = 1;
    for i in 2..=n as u64 {
        acc = acc
            .checked_mul(i)
            .ok_or_else(|| OrbitalError::Numerical("factorial overflow".into()))?;
    }

    Ok(acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn laguerre_basic_values() {
        assert_relative_eq!(associated_laguerre(0, 1, 2.0).unwrap(), 1.0);
        assert_relative_eq!(associated_laguerre(1, 1, 2.0).unwrap(), 0.0);
    }

    #[test]
    fn one_s_at_origin_is_finite() {
        let cfg = OrbitalConfig::ground_state(1.0, 1.0).unwrap();
        let r = RadialWavefunction::new(cfg);
        let val = r.value(0.0).unwrap();

        assert!(val.is_finite());
        assert!(val > 0.0);
    }

    #[test]
    fn three_d_z2_radial_node_at_origin() {
        let cfg = OrbitalConfig::d_z2(1.0, 1.0).unwrap();
        let r = RadialWavefunction::new(cfg);
        let val = r.value(0.0).unwrap();

        assert_relative_eq!(val, 0.0, epsilon = 1e-14);
    }
}
