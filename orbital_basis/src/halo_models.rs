// horb/orbital_basis/src/halo_models.rs

use crate::{OrbitalError, G_KPC};
use statrs::function::gamma::{gamma_ur, ln_gamma};
use std::f64::consts::PI;

/// Default for reproducible fixed-profile comparisons, not a fitted value.
pub const DEFAULT_EINASTO_ALPHA: f64 = 0.18;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClassicalHaloKind {
    PseudoIsothermal,
    Nfw,
    Burkert,
    Einasto,
}

impl ClassicalHaloKind {
    pub fn parse(s: &str) -> Result<Self, OrbitalError> {
        match s {
            "piso" | "pseudo_isothermal" | "pseudo-isothermal" => Ok(Self::PseudoIsothermal),
            "nfw" | "NFW" => Ok(Self::Nfw),
            "burkert" => Ok(Self::Burkert),
            "einasto" | "Einasto" => Ok(Self::Einasto),
            _ => Err(OrbitalError::InvalidParameter {
                name: "classical_halo_kind",
                reason: format!("unknown halo kind '{s}'; use piso, nfw, burkert, or einasto"),
            }),
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            Self::PseudoIsothermal => "piso",
            Self::Nfw => "nfw",
            Self::Burkert => "burkert",
            Self::Einasto => "einasto",
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
///   core radius for Burkert, slope-minus-two radius r_-2 for Einasto.
///
/// density_norm:
///   rho0 for pseudo-isothermal and Burkert,
///   rho_s for NFW, rho_-2 for Einasto.
///
/// The constructor `from_mass_at_radius` solves density_norm so that:
///
///   M(<r_ref_kpc) = m_ref_msun
#[derive(Debug, Clone, Copy)]
pub struct ClassicalHalo {
    pub kind: ClassicalHaloKind,
    pub scale_kpc: f64,
    pub density_norm_msun_per_kpc3: f64,
    /// Einasto curvature alpha = 1/n. Ignored by the other halo families.
    /// Use `einasto_from_mass_at_radius` to set a nondefault value.
    pub einasto_alpha: f64,
}

impl ClassicalHalo {
    pub fn from_mass_at_radius(
        kind: ClassicalHaloKind,
        scale_kpc: f64,
        m_ref_msun: f64,
        r_ref_kpc: f64,
    ) -> Result<Self, OrbitalError> {
        Self::with_shape(
            kind,
            scale_kpc,
            m_ref_msun,
            r_ref_kpc,
            DEFAULT_EINASTO_ALPHA,
        )
    }

    /// rho(r) = rho_-2 exp[-2/alpha ((r/r_-2)^alpha - 1)].
    /// The supplied mass is enclosed within r_ref, not the mass at infinity.
    /// See Retana-Montenegro et al. (2012), A&A 540, A70, Eqs. 4, 8, 12-13:
    /// <https://doi.org/10.1051/0004-6361/201118543>.
    pub fn einasto_from_mass_at_radius(
        scale_kpc: f64,
        alpha: f64,
        m_ref_msun: f64,
        r_ref_kpc: f64,
    ) -> Result<Self, OrbitalError> {
        Self::with_shape(
            ClassicalHaloKind::Einasto,
            scale_kpc,
            m_ref_msun,
            r_ref_kpc,
            alpha,
        )
    }

    fn with_shape(
        kind: ClassicalHaloKind,
        scale_kpc: f64,
        m_ref_msun: f64,
        r_ref_kpc: f64,
        alpha: f64,
    ) -> Result<Self, OrbitalError> {
        if !scale_kpc.is_finite() || scale_kpc <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "scale_kpc",
                reason: "must be finite and positive".into(),
            });
        }

        if !m_ref_msun.is_finite() || m_ref_msun <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "m_ref_msun",
                reason: "must be finite and positive".into(),
            });
        }

        if !r_ref_kpc.is_finite() || r_ref_kpc <= 0.0 {
            return Err(OrbitalError::InvalidParameter {
                name: "r_ref_kpc",
                reason: "must be finite and positive".into(),
            });
        }

        if !alpha.is_finite() || alpha <= 0.0 || !(3.0 / alpha).is_finite() {
            return Err(OrbitalError::InvalidParameter {
                name: "einasto_alpha",
                reason: "must be finite, positive and representable in f64".into(),
            });
        }

        let density_norm_msun_per_kpc3 = if kind == ClassicalHaloKind::Einasto {
            (m_ref_msun.ln() - einasto_log_unit_mass(scale_kpc, alpha, r_ref_kpc)).exp()
        } else {
            m_ref_msun / Self::unit_density_mass(kind, scale_kpc, r_ref_kpc)
        };
        if !density_norm_msun_per_kpc3.is_finite() || density_norm_msun_per_kpc3 <= 0.0 {
            return Err(OrbitalError::Numerical(
                "halo normalization is not positive and finite; check parameter magnitudes".into(),
            ));
        }

        Ok(Self {
            kind,
            scale_kpc,
            density_norm_msun_per_kpc3,
            einasto_alpha: alpha,
        })
    }

    pub fn density(&self, r_kpc: f64) -> f64 {
        if r_kpc.is_nan() || r_kpc < 0.0 {
            return f64::NAN;
        }

        if self.kind == ClassicalHaloKind::Einasto {
            // expm1 preserves the slope as alpha becomes small or r approaches r_-2.
            let log_x = r_kpc.ln() - self.scale_kpc.ln();
            return (self.density_norm_msun_per_kpc3.ln()
                - (2.0 / self.einasto_alpha) * (self.einasto_alpha * log_x).exp_m1())
            .exp();
        }

        let r = r_kpc.max(1.0e-12);
        let x = r / self.scale_kpc;
        let rho = self.density_norm_msun_per_kpc3;

        match self.kind {
            ClassicalHaloKind::PseudoIsothermal => rho / (1.0 + x * x),
            ClassicalHaloKind::Nfw => rho / (x * (1.0 + x).powi(2)),
            ClassicalHaloKind::Burkert => rho / ((1.0 + x) * (1.0 + x * x)),
            ClassicalHaloKind::Einasto => unreachable!(),
        }
    }

    pub fn enclosed_mass(&self, r_kpc: f64) -> f64 {
        if r_kpc <= 0.0 {
            return 0.0;
        }

        if self.kind == ClassicalHaloKind::Einasto {
            return (self.density_norm_msun_per_kpc3.ln()
                + einasto_log_unit_mass(self.scale_kpc, self.einasto_alpha, r_kpc))
            .exp();
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
            ClassicalHaloKind::Einasto => unreachable!("Einasto uses log-space mass integration"),
        }
    }
}

/// ln M(<r) with rho_-2 = 1. The small-x branch combines the lower-gamma
/// series with the prefactor analytically, avoiding both cancellation and
/// underflow of the regularized gamma at small alpha/radius.
fn einasto_log_unit_mass(scale: f64, alpha: f64, r: f64) -> f64 {
    if r.is_nan() {
        return f64::NAN;
    }
    if r <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let s = 3.0 / alpha;
    let b = 2.0 / alpha;
    let log_x = r.ln() - scale.ln();
    let t = (b.ln() + alpha * log_x).exp();
    if t < s + 1.0 {
        // gamma(s,t) = t^s exp(-t) / s *
        //              [1 + t/(s+1) + t^2/((s+1)(s+2)) + ...].
        let mut term = 1.0;
        let mut sum = 1.0;
        for k in 1..=100_000 {
            term *= t / (s + k as f64);
            sum += term;
            if term <= 2.0e-15 * sum {
                return (4.0 * PI / 3.0).ln() + 3.0 * r.ln() - b * (alpha * log_x).exp_m1()
                    + sum.ln();
            }
        }
        return f64::NAN; // Explicit failure, not a silently unconverged value.
    }
    let log_p = if t.is_infinite() {
        0.0
    } else {
        (-gamma_ur(s, t)).ln_1p()
    };
    (4.0 * PI).ln() + 3.0 * scale.ln() + b - s * b.ln() - alpha.ln() + ln_gamma(s) + log_p
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
            ClassicalHaloKind::Einasto,
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

    #[test]
    fn einasto_matches_independent_high_precision_reference_values() {
        // mpmath 1.3.0, 70 decimal digits, lower incomplete gamma expression.
        // r_-2=15 kpc, M(<80 kpc)=2.5e11 Msun. alpha=0.005 has Gamma(600),
        // which would overflow if the gamma and prefactor were multiplied directly.
        for (alpha, rho, m8, m200) in [
            (
                0.005,
                1_113_253.310_852_480_7,
                24_736_964_962.348_423,
                618_525_689_288.043,
            ),
            (
                0.05,
                1_184_622.972_504_932_5,
                23_067_113_744.148_785,
                561_973_379_860.754_3,
            ),
            (
                0.18,
                1_401_820.424_937_787,
                21_128_339_371.231_983,
                427_302_034_831.128_36,
            ),
            (
                0.5,
                2_051_885.839_991_037_8,
                21_208_828_227.845_93,
                277_377_832_341.414_3,
            ),
            (
                1.0,
                3_196_111.626_327_174_4,
                23_278_091_419.341_503,
                250_400_131_240.351,
            ),
            (
                2.0,
                4_893_808.157_716_53,
                24_120_968_949.007_637,
                250_000_000_000.678_53,
            ),
            (
                5.0,
                7_655_856.892_768_549,
                24_336_882_357.833_992,
                250_000_000_000.0,
            ),
        ] {
            let h = ClassicalHalo::einasto_from_mass_at_radius(15.0, alpha, 2.5e11, 80.0).unwrap();
            for (actual, expected) in [
                (h.density(15.0), rho),
                (h.enclosed_mass(8.0), m8),
                (h.enclosed_mass(200.0), m200),
            ] {
                assert!(
                    (actual / expected - 1.0).abs() < 2e-11,
                    "alpha={alpha}: {actual} vs {expected}"
                );
            }
        }
    }

    #[test]
    fn einasto_mass_matches_density_quadrature() {
        // Simpson integration in u, with r=R*u^4 to resolve the central density.
        for alpha in [0.05, 0.18, 0.5, 1.0, 2.0] {
            let h = ClassicalHalo::einasto_from_mass_at_radius(15.0, alpha, 2.5e11, 80.0).unwrap();
            for radius in [0.05_f64, 8.0, 80.0, 200.0] {
                let n = 16_000;
                let mut sum = 0.0;
                for i in 1..=n {
                    let u = i as f64 / n as f64;
                    let weight = if i == n {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    sum += weight
                        * 16.0
                        * PI
                        * radius.powi(3)
                        * u.powi(11)
                        * h.density(radius * u.powi(4));
                }
                let numerical = sum / (3.0 * n as f64);
                assert!(
                    (h.enclosed_mass(radius) / numerical - 1.0).abs() < 2e-8,
                    "alpha={alpha}, R={radius}"
                );
            }
        }
    }

    #[test]
    fn einasto_exponential_special_case_and_small_radius_limit() {
        let h = ClassicalHalo::einasto_from_mass_at_radius(10.0, 1.0, 1e11, 80.0).unwrap();
        for r in [1.0_f64, 10.0, 50.0, 200.0] {
            let x = r / 10.0;
            let expected = PI
                * h.density(0.0)
                * 10.0_f64.powi(3)
                * (1.0 - (-2.0 * x).exp() * (1.0 + 2.0 * x + 2.0 * x * x));
            assert!((h.enclosed_mass(r) / expected - 1.0).abs() < 2e-12);
        }
        let r = 1e-12_f64;
        let core = 4.0 * PI / 3.0 * h.density(0.0) * r.powi(3);
        assert!((h.enclosed_mass(r) / core - 1.0).abs() < 1e-11);
        assert_eq!(h.enclosed_mass(0.0), 0.0);
        assert_eq!(h.circular_velocity(0.0), 0.0);
    }

    #[test]
    fn einasto_scale_is_slope_minus_two_and_mass_derivative_matches_density() {
        let h = ClassicalHalo::einasto_from_mass_at_radius(15.0, 0.18, 2.5e11, 80.0).unwrap();
        let eps = 1e-5_f64;
        let slope =
            (h.density(15.0 * eps.exp()).ln() - h.density(15.0 * (-eps).exp()).ln()) / (2.0 * eps);
        assert!((slope + 2.0).abs() < 1e-8);
        for r in [0.01_f64, 1.0, 15.0, 80.0, 300.0] {
            let dr = r * 1e-5;
            let derivative = (h.enclosed_mass(r + dr) - h.enclosed_mass(r - dr)) / (2.0 * dr);
            assert!((derivative / (4.0 * PI * r * r * h.density(r)) - 1.0).abs() < 2e-7);
        }
    }

    #[test]
    fn einasto_has_finite_total_mass_and_correct_mass_scaling() {
        let h = ClassicalHalo::einasto_from_mass_at_radius(15.0, 0.18, 2.5e11, 80.0).unwrap();
        let twice = ClassicalHalo::einasto_from_mass_at_radius(15.0, 0.18, 5e11, 80.0).unwrap();
        let total = h.enclosed_mass(f64::INFINITY);
        assert!(total.is_finite() && total > 2.5e11);
        assert!((h.enclosed_mass(1e9) / total - 1.0).abs() < 1e-10);
        assert!((twice.enclosed_mass(8.0) / h.enclosed_mass(8.0) - 2.0).abs() < 1e-12);
        assert!(
            (twice.circular_velocity(8.0) / h.circular_velocity(8.0) - 2.0_f64.sqrt()).abs()
                < 1e-12
        );
        assert_eq!(h.density(f64::INFINITY), 0.0);
        assert!(h.density(f64::NAN).is_nan());
        assert!(h.enclosed_mass(f64::NAN).is_nan());
    }

    #[test]
    fn einasto_rejects_invalid_parameters_and_parses_its_name() {
        assert_eq!(
            ClassicalHaloKind::parse("einasto").unwrap().label(),
            "einasto"
        );
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(ClassicalHalo::einasto_from_mass_at_radius(15.0, bad, 2.5e11, 80.0).is_err());
            assert!(ClassicalHalo::einasto_from_mass_at_radius(bad, 0.18, 2.5e11, 80.0).is_err());
            assert!(ClassicalHalo::einasto_from_mass_at_radius(15.0, 0.18, bad, 80.0).is_err());
            assert!(ClassicalHalo::einasto_from_mass_at_radius(15.0, 0.18, 2.5e11, bad).is_err());
        }
    }
}
