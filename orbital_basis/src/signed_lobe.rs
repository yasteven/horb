//! Real wavefunction used as a phenomenological signed mass density.
//!
//! rho = M_ref psi / integral(|psi| dV), NOT a quantum probability density.
//! Normalize by absolute mass, since every l>0 state has zero net mass.
//! The real phase convention is the one in `angular` and `radial`.
use crate::{
    signed_gravity::{gauss_legendre, MassCell},
    OrbitalConfig, OrbitalError, RadialWavefunction, RealSphericalHarmonic, G_KPC,
};
use std::f64::consts::PI;

#[derive(Clone)]
pub struct SignedLobe {
    cfg: OrbitalConfig,
    radial: RadialWavefunction,
    angular: RealSphericalHarmonic,
    nodes: Vec<(f64, f64)>,
    breaks: Vec<f64>,
    norm: f64,
    pub net_mass_fraction: f64,
}

impl SignedLobe {
    pub fn new(cfg: OrbitalConfig, order: usize) -> Result<Self, OrbitalError> {
        crate::signed_gravity::state_config(
            &format!("{}:{}:{}", cfg.n, cfg.l, cfg.m),
            cfg.a0_star,
            cfg.dm_mass,
        )?;
        let mut unit = cfg.clone();
        unit.a0_star = 1.0;
        unit.dm_mass = 1.0;
        let radial = RadialWavefunction::new(unit);
        let nodes = gauss_legendre(order)?;
        let limit = (16.0 * (cfg.n as f64).powi(2)).max(120.0 * cfg.n as f64);
        // Split at radial nodes before integrating |R|. A single quadrature
        // across the cusps would give order-dependent mass normalization.
        let mut breaks = vec![0.0];
        let mut previous = limit / 4096.0;
        let mut value = radial.value(previous)?;
        for i in 2..=4096 {
            let next = limit * i as f64 / 4096.0;
            let next_value = radial.value(next)?;
            if value.signum() != next_value.signum() && value != 0.0 && next_value != 0.0 {
                let (mut lo, mut hi) = (previous, next);
                for _ in 0..60 {
                    let mid = (lo + hi) / 2.0;
                    if radial.value(mid)?.signum() == value.signum() {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                breaks.push((lo + hi) / 2.0);
            }
            previous = next;
            value = next_value;
        }
        if breaks.len() != (cfg.n - cfg.l) as usize {
            return Err(OrbitalError::Numerical(
                "signed-lobe radial node isolation failed".into(),
            ));
        }
        breaks.push(limit);
        let angular_abs = match (cfg.l, cfg.m) {
            (0, _) => (4.0 * PI).sqrt(),
            (1, _) => (3.0 * PI).sqrt(),
            (2, 0) => 4.0 * (5.0 * PI).sqrt() / (3.0 * 3.0_f64.sqrt()),
            (2, _) => 4.0 * (15.0 / PI).sqrt() / 3.0,
            _ => unreachable!(),
        };
        let angular = RealSphericalHarmonic::new(cfg.l, cfg.m)?;
        let mut field = Self {
            cfg,
            radial,
            angular,
            nodes,
            breaks,
            norm: 1.0,
            net_mass_fraction: 0.0,
        };
        field.norm = angular_abs * field.integrate(0.0, limit, 2, true)?;
        if field.cfg.l == 0 {
            field.net_mass_fraction =
                (4.0 * PI).sqrt() * field.integrate(0.0, limit, 2, false)? / field.norm;
        }
        Ok(field)
    }

    fn integrate(&self, lo: f64, hi: f64, power: i32, absolute: bool) -> Result<f64, OrbitalError> {
        let mut sum = 0.0;
        for interval in self.breaks.windows(2) {
            let left = lo.max(interval[0]);
            let right = hi.min(interval[1]);
            if right <= left {
                continue;
            }
            let half = (right - left) / 2.0;
            let mid = (right + left) / 2.0;
            for (node, weight) in &self.nodes {
                let r = mid + half * node;
                let mut value = self.radial.value(r)?;
                if absolute {
                    value = value.abs();
                }
                sum += half * weight * value * r.powi(power);
            }
        }
        Ok(sum)
    }

    pub fn rho(&self, r: f64, theta: f64, phi: f64) -> Result<f64, OrbitalError> {
        Ok(self.cfg.dm_mass / self.cfg.a0_star.powi(3) / self.norm
            * self.radial.value(r / self.cfg.a0_star)?
            * self.angular.value(theta, phi)?)
    }

    /// Azimuthally averaged R*dPhi/dR in z=0, with full exterior multipoles.
    /// Odd l and m!=0 have exactly zero response in this observable.
    pub fn disk_plane_v2(&self, r: f64) -> Result<f64, OrbitalError> {
        if !r.is_finite() || r < 0.0 {
            return Err(OrbitalError::Numerical(
                "radius must be finite and nonnegative".into(),
            ));
        }
        if r == 0.0 || self.cfg.m != 0 || self.cfg.l % 2 == 1 {
            return Ok(0.0);
        }
        let x = r / self.cfg.a0_star;
        let ell = self.cfg.l as i32;
        let inside = self.integrate(0.0, x, ell + 2, false)?;
        let outside = if ell == 0 {
            0.0
        } else {
            self.integrate(x, *self.breaks.last().unwrap(), 1 - ell, false)?
        };
        let equator = if ell == 0 { 1.0 } else { -0.5 };
        let angular = (4.0 * PI / (2 * ell + 1) as f64).sqrt() * equator / self.norm;
        Ok(G_KPC * self.cfg.dm_mass / self.cfg.a0_star
            * angular
            * ((ell + 1) as f64 * inside / x.powi(ell + 1) - ell as f64 * x.powi(ell) * outside))
    }

    /// Direct-force quadrature for the axisymmetric model (CPU/CUDA audit).
    pub fn mass_cells(&self, nmu: usize, nphi: usize) -> Result<Vec<MassCell>, OrbitalError> {
        if self.cfg.m != 0 || !(8..=512).contains(&nphi) {
            return Err(OrbitalError::Numerical(
                "signed-lobe force cells require m=0 and n_phi=8..512".into(),
            ));
        }
        let angular = gauss_legendre(nmu)?;
        if (self.breaks.len() - 1) * self.nodes.len() * nmu * nphi > 8_000_000 {
            return Err(OrbitalError::Numerical(
                "at most eight million mass cells".into(),
            ));
        }
        let mut cells = Vec::new();
        for interval in self.breaks.windows(2) {
            let half = (interval[1] - interval[0]) * self.cfg.a0_star / 2.0;
            let mid = (interval[1] + interval[0]) * self.cfg.a0_star / 2.0;
            for (node, wr) in &self.nodes {
                let r = mid + half * node;
                for &(mu, wm) in &angular {
                    let mass = self.rho(r, mu.acos(), 0.0)? * r * r * half * wr * wm * 2.0 * PI
                        / nphi as f64;
                    for j in 0..nphi {
                        let phi = 2.0 * PI * (j as f64 + 0.5) / nphi as f64;
                        let s = (1.0 - mu * mu).sqrt();
                        cells.push(MassCell {
                            x: r * s * phi.cos(),
                            y: r * s * phi.sin(),
                            z: r * mu,
                            m: mass,
                        });
                    }
                }
            }
        }
        Ok(cells)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signed_gravity::{direct_v2, state_config};
    #[test]
    fn one_s_matches_exponential_mass_and_scaling() {
        let field = SignedLobe::new(state_config("1s", 2.0, 1e11).unwrap(), 128).unwrap();
        for r in [0.1, 2.0, 10.0, 80.0] {
            let x: f64 = r / 2.0;
            let mass = 1e11 * (1.0 - (-x).exp() * (1.0 + x + x * x / 2.0));
            assert!((field.disk_plane_v2(r).unwrap() / (G_KPC * mass / r) - 1.0).abs() < 1e-9);
        }
        assert!((field.net_mass_fraction - 1.0).abs() < 1e-12);
    }
    #[test]
    fn signs_nodes_and_zero_net_mass_are_retained() {
        let field = SignedLobe::new(state_config("3d_z2", 1.0, 1e11).unwrap(), 128).unwrap();
        assert!(field.rho(3.0, 0.0, 0.0).unwrap() > 0.0);
        assert!(field.rho(3.0, PI / 2.0, 0.0).unwrap() < 0.0);
        assert_eq!(field.net_mass_fraction, 0.0);
        let cells = field.mass_cells(32, 64).unwrap();
        assert!(cells.iter().map(|c| c.m).sum::<f64>().abs() < 0.01);
        for r in [80.0, 160.0] {
            let exact = field.disk_plane_v2(r).unwrap();
            assert!((direct_v2(&cells, r, 0.0) / exact - 1.0).abs() < 2e-4);
        }
        let s = SignedLobe::new(state_config("2s", 1.0, 1e11).unwrap(), 128).unwrap();
        assert!(s.rho(1.0, 0.0, 0.0).unwrap() > 0.0);
        assert!(s.rho(3.0, 0.0, 0.0).unwrap() < 0.0);
        assert!(s.net_mass_fraction < 0.0);
        let p = SignedLobe::new(state_config("2p_z", 1.0, 1e11).unwrap(), 128).unwrap();
        assert_eq!(p.disk_plane_v2(2.0).unwrap(), 0.0);
    }
    #[test]
    fn higher_states_converge_with_absolute_normalization() {
        for n in 1..=10 {
            for l in 0..n.min(3) {
                let cfg = state_config(&format!("{n}:{l}:0"), 1.0, 1.0).unwrap();
                let a = SignedLobe::new(cfg.clone(), 128).unwrap();
                let b = SignedLobe::new(cfg, 256).unwrap();
                assert!((a.norm / b.norm - 1.0).abs() < 1e-10);
                for r in [0.1, 1.0, 10.0, 100.0, 1000.0] {
                    assert!(
                        (a.disk_plane_v2(r).unwrap() - b.disk_plane_v2(r).unwrap()).abs() < 1e-14
                    );
                }
            }
        }
    }
}
