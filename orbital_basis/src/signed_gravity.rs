//! Linear disk-plane gravity of hydrogenic density basis functions.
//! The CPU multipole result is for the azimuthal mean, with no softening.
//! For l <= 2, |Y_lm|^2 needs only L=0,2,4; this is not a spherical mass
//! approximation. Nonaxisymmetric orbitals retain their azimuthal mean only.
use crate::{OrbitalConfig, OrbitalError, RadialWavefunction, G_KPC};
use std::f64::consts::PI;

fn invalid(name: &'static str, reason: &str) -> OrbitalError {
    OrbitalError::InvalidParameter {
        name,
        reason: reason.into(),
    }
}

pub fn state_config(name: &str, scale: f64, mass: f64) -> Result<OrbitalConfig, OrbitalError> {
    let q = match name {
        "1s" => (1, 0, 0),
        "2s" => (2, 0, 0),
        "3s" => (3, 0, 0),
        "2p_z" | "2pz" => (2, 1, 0),
        "2p_x" => (2, 1, 1),
        "2p_y" => (2, 1, -1),
        "3p_z" | "3pz" => (3, 1, 0),
        "3p_x" => (3, 1, 1),
        "3p_y" => (3, 1, -1),
        "3d_z2" => (3, 2, 0),
        "3d_xz" => (3, 2, 1),
        "3d_yz" => (3, 2, -1),
        "3d_x2_y2" => (3, 2, 2),
        "3d_xy" => (3, 2, -2),
        _ => {
            let parts: Vec<_> = name.split(':').collect();
            if parts.len() != 3 {
                return Err(invalid("state", "use a named orbital or n:l:m"));
            }
            (
                parts[0]
                    .parse()
                    .map_err(|_| invalid("state", "invalid n"))?,
                parts[1]
                    .parse()
                    .map_err(|_| invalid("state", "invalid l"))?,
                parts[2]
                    .parse()
                    .map_err(|_| invalid("state", "invalid m"))?,
            )
        }
    };
    let cfg = OrbitalConfig::new(q.0, q.1, q.2, scale, mass)?;
    validate(&cfg)?;
    Ok(cfg)
}

fn validate(cfg: &OrbitalConfig) -> Result<(), OrbitalError> {
    OrbitalConfig::new(cfg.n, cfg.l, cfg.m, cfg.a0_star, cfg.dm_mass)?;
    if cfg.l > 2 || cfg.n > 10 {
        return Err(invalid("state", "signed basis supports n <= 10 and l <= 2"));
    }
    Ok(())
}

/// Gauss-Legendre nodes/weights on [-1,1].
pub fn gauss_legendre(n: usize) -> Result<Vec<(f64, f64)>, OrbitalError> {
    if !(8..=512).contains(&n) {
        return Err(invalid("quadrature order", "must be 8..512"));
    }
    let mut nodes = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 0.0;
        let mut converged = false;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let pk = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = pk;
            }
            derivative = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let step = p1 / derivative;
            x -= step;
            if step.abs() < 2e-15 {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(OrbitalError::Numerical(
                "Gauss-Legendre root did not converge".into(),
            ));
        }
        nodes.push((x, 2.0 / ((1.0 - x * x) * derivative * derivative)));
    }
    nodes.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(nodes)
}

fn legendre(l: u32, mu: f64) -> f64 {
    match l {
        0 => 1.0,
        2 => (3.0 * mu * mu - 1.0) / 2.0,
        4 => (35.0 * mu.powi(4) - 30.0 * mu * mu + 3.0) / 8.0,
        _ => unreachable!(),
    }
}

/// Azimuthal average of the normalized real harmonic's squared amplitude.
pub fn angular_mean(l: u32, m: i32, mu: f64) -> f64 {
    let s2 = 1.0 - mu * mu;
    match (l, m.unsigned_abs()) {
        (0, 0) => 1.0 / (4.0 * PI),
        (1, 0) => 3.0 / (4.0 * PI) * mu * mu,
        (1, 1) => 3.0 / (8.0 * PI) * s2,
        (2, 0) => 5.0 / (16.0 * PI) * (3.0 * mu * mu - 1.0).powi(2),
        (2, 1) => 15.0 / (8.0 * PI) * s2 * mu * mu,
        (2, 2) => 15.0 / (32.0 * PI) * s2 * s2,
        _ => f64::NAN,
    }
}

/// Dimensionless radial cutoff with an exponentially negligible tail for n<=10.
fn cutoff(cfg: &OrbitalConfig) -> f64 {
    (8.0 * (cfg.n as f64).powi(2)).max(60.0 * cfg.n as f64)
}

pub fn disk_plane_v2(cfg: &OrbitalConfig, radius: f64, order: usize) -> Result<f64, OrbitalError> {
    validate(cfg)?;
    if !radius.is_finite() || radius < 0.0 {
        return Err(invalid("radius", "must be finite and nonnegative"));
    }
    let nodes = gauss_legendre(order)?;
    if radius == 0.0 {
        return Ok(0.0);
    }
    // Work in r/a0 to avoid dimensional overflow and resolve compact states.
    let mut unit_cfg = cfg.clone();
    unit_cfg.a0_star = 1.0;
    unit_cfg.dm_mass = 1.0;
    let radial = RadialWavefunction::new(unit_cfg);
    let x = radius / cfg.a0_star;
    let limit = cutoff(cfg);
    let integrate = |lo: f64, hi: f64, power: i32| -> Result<f64, OrbitalError> {
        if hi <= lo {
            return Ok(0.0);
        }
        let mid = (hi + lo) * 0.5;
        let half = (hi - lo) * 0.5;
        nodes.iter().try_fold(0.0, |sum, (node, weight)| {
            let u = mid + half * node;
            Ok(sum + half * weight * radial.radial_density_factor(u)? * u.powi(power))
        })
    };
    let angular_nodes = gauss_legendre(8)?;
    let mut k = 0.0;
    for ell in (0..=2 * cfg.l).step_by(2) {
        let a: f64 = (2 * ell + 1) as f64 / 2.0
            * angular_nodes
                .iter()
                .map(|(mu, w)| w * angular_mean(cfg.l, cfg.m, *mu) * legendre(ell, *mu))
                .sum::<f64>();
        let inside = integrate(0.0, x.min(limit), ell as i32 + 2)?;
        let outside = if ell == 0 {
            0.0
        } else {
            integrate(x, limit, 1 - ell as i32)?
        };
        k += 4.0 * PI / (2 * ell + 1) as f64
            * a
            * legendre(ell, 0.0)
            * ((ell + 1) as f64 * inside / x.powi(ell as i32 + 1)
                - ell as f64 * x.powi(ell as i32) * outside);
    }
    let result = G_KPC * cfg.dm_mass / cfg.a0_star * k;
    if !result.is_finite() {
        return Err(OrbitalError::Numerical("nonfinite disk-plane force".into()));
    }
    Ok(result)
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct MassCell {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub m: f64,
}

/// Independent quadrature cells for CUDA/direct-force verification. No mass
/// renormalization: weights retain the full-space analytic orbital mass.
pub fn mass_cells(
    cfg: &OrbitalConfig,
    nr: usize,
    nmu: usize,
    nphi: usize,
) -> Result<Vec<MassCell>, OrbitalError> {
    validate(cfg)?;
    if !(8..=512).contains(&nphi) {
        return Err(invalid("n_phi", "must be 8..512"));
    }
    let radial_nodes = gauss_legendre(nr)?;
    let angular_nodes = gauss_legendre(nmu)?;
    if nr * nmu * nphi > 8_000_000 {
        return Err(invalid("quadrature", "at most eight million cells"));
    }
    let radial = RadialWavefunction::new(cfg.clone());
    let limit = cutoff(cfg) * cfg.a0_star;
    let mut cells = Vec::with_capacity(nr * nmu * nphi);
    for (u, wr) in radial_nodes {
        let r = (u + 1.0) * limit * 0.5;
        let radial_mass = cfg.dm_mass * radial.radial_density_factor(r)? * r * r * wr * limit * 0.5;
        for &(mu, wm) in &angular_nodes {
            let m = radial_mass * angular_mean(cfg.l, cfg.m, mu) * wm * 2.0 * PI / nphi as f64;
            for j in 0..nphi {
                let phi = 2.0 * PI * (j as f64 + 0.5) / nphi as f64;
                let s = (1.0 - mu * mu).sqrt();
                cells.push(MassCell {
                    x: r * s * phi.cos(),
                    y: r * s * phi.sin(),
                    z: r * mu,
                    m,
                });
            }
        }
    }
    Ok(cells)
}

/// Signed radial force at (R,0,0), multiplied by R. Never clamp or sqrt.
pub fn direct_v2(cells: &[MassCell], r: f64, softening: f64) -> f64 {
    G_KPC
        * r
        * cells
            .iter()
            .map(|c| {
                let dx = r - c.x;
                let d2 = dx * dx + c.y * c.y + c.z * c.z + softening * softening;
                c.m * dx / (d2 * d2.sqrt())
            })
            .sum::<f64>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multipole_matches_analytic_one_s_and_has_correct_scaling() {
        let cfg = state_config("1s", 2.0, 1e11).unwrap();
        for r in [0.01, 0.5, 2.0, 10.0, 80.0] {
            let x: f64 = r / 2.0;
            let mass = 1e11 * (1.0 - (-2.0 * x).exp() * (1.0 + 2.0 * x + 2.0 * x * x));
            let k = disk_plane_v2(&cfg, r, 128).unwrap();
            assert!((k / (G_KPC * mass / r) - 1.0).abs() < 1e-8, "r={r}");
            let scaled = state_config("1s", 4.0, 2e11).unwrap();
            assert!((disk_plane_v2(&scaled, 2.0 * r, 128).unwrap() / k - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn angular_means_and_cells_integrate_to_reference_mass() {
        for name in [
            "1s", "2s", "2p_z", "2p_x", "3s", "3p_z", "3d_z2", "3d_xy", "3d_xz",
        ] {
            let cfg = state_config(name, 1.0, 1e11).unwrap();
            let cells = mass_cells(&cfg, 128, 16, 16).unwrap();
            assert!(
                (cells.iter().map(|c| c.m).sum::<f64>() / 1e11 - 1.0).abs() < 1e-9,
                "{name}"
            );
        }
    }

    #[test]
    fn multipole_agrees_with_independent_direct_force_outside_halo() {
        for name in ["2p_z", "3d_z2", "3d_xy", "3d_xz"] {
            let cfg = state_config(name, 1.0, 1e11).unwrap();
            let cells = mass_cells(&cfg, 128, 32, 64).unwrap();
            for r in [40.0, 80.0] {
                let direct = direct_v2(&cells, r, 0.0);
                assert!(
                    (disk_plane_v2(&cfg, r, 128).unwrap() / direct - 1.0).abs() < 5e-4,
                    "{name} R={r}"
                );
            }
        }
    }

    #[test]
    fn signed_cells_preserve_outward_force() {
        let cells = [
            MassCell {
                x: 10.0,
                y: 0.0,
                z: 0.0,
                m: 1e10,
            },
            MassCell {
                x: -10.0,
                y: 0.0,
                z: 0.0,
                m: 1e10,
            },
        ];
        assert!(direct_v2(&cells, 1.0, 0.1) < 0.0);
        let negative: Vec<_> = cells.iter().map(|c| MassCell { m: -c.m, ..*c }).collect();
        assert_eq!(direct_v2(&negative, 1.0, 0.1), -direct_v2(&cells, 1.0, 0.1));
    }

    #[test]
    fn nonspherical_inner_force_matches_high_precision_multipoles() {
        // Independent mpmath incomplete-gamma evaluation (60 digits) of the
        // 3d_z2 L=0,2,4 radial moments. This checks both interior/exterior terms.
        let cfg = state_config("3d_z2", 1.0, 1e11).unwrap();
        for (r, expected) in [
            (0.05, 0.7583256035664207),
            (1.0, 275.8974661279248),
            (5.0, 5600.430683727943),
            (10.0, 20447.645268698514),
            (25.0, 16259.64979783057),
        ] {
            assert!((disk_plane_v2(&cfg, r, 128).unwrap() / expected - 1.0).abs() < 1e-10);
        }
        let torus = state_config("3d_xy", 1.0, 1e11).unwrap();
        assert!(disk_plane_v2(&torus, 0.5, 128).unwrap() < 0.0);
        assert!(disk_plane_v2(&torus, 80.0, 128).unwrap() > 0.0);
    }
}
