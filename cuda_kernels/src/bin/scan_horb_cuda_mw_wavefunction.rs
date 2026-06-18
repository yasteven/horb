use anyhow::{bail, Result};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct OrbitalTerm {
    pub state_id: i32,
    pub coeff: f64,
}

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct EulerRotation {
    pub r00: f64,
    pub r01: f64,
    pub r02: f64,
    pub r10: f64,
    pub r11: f64,
    pub r12: f64,
    pub r20: f64,
    pub r21: f64,
    pub r22: f64,
}

#[link(name = "horb_cuda", kind = "dylib")]
unsafe extern "C" {
    fn compute_superposition_disk_curve_cuda(
        terms: *const OrbitalTerm,
        term_count: i32,
        a0: f64,
        total_mass_msun: f64,
        n_side: i32,
        extent_kpc: f64,
        softening_kpc: f64,
        radii: *const f64,
        v_out: *mut f64,
        radius_count: usize,
        rot: EulerRotation,
    );

    fn make_identity_rotation_cuda() -> EulerRotation;
}

#[derive(Debug, Clone)]
struct RcRow {
    r_kpc: f64,
    vobs_kms: f64,
    ev_kms: f64,
}

#[derive(Debug, Clone)]
struct BaryonComponent {
    kind: String,
    mass_msun: f64,
    scale_kpc: f64,
}

const G_KPC: f64 = 4.3009e-6;

fn state_id(name: &str) -> Result<i32> {
    match name {
        "1s" => Ok(0),
        "2p_z" => Ok(10),
        "2p_x" => Ok(11),
        "2p_y" => Ok(12),
        "3d_z2" => Ok(20),
        "3d_x2_y2" => Ok(21),
        "3d_xy" => Ok(22),
        "3d_xz" => Ok(23),
        "3d_yz" => Ok(24),
        _ => bail!("unknown state: {name}"),
    }
}

fn parse_list_f64(value: &str) -> Result<Vec<f64>> {
    value.split(',').map(|v| Ok(v.trim().parse::<f64>()?)).collect()
}

fn parse_list_string(value: &str) -> Vec<String> {
    value
        .split(',')
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .collect()
}

fn parse_standard_rc(path: &str, r_min: f64, r_max: f64) -> Result<Vec<RcRow>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    let mut rows = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        if line_no == 0 && line.starts_with("R_kpc") {
            continue;
        }

        let parts: Vec<&str> = line.split(',').map(|v| v.trim()).collect();

        if parts.len() < 3 {
            bail!("bad RC row {} in {}", line_no + 1, path);
        }

        let r_kpc: f64 = parts[0].parse()?;
        let vobs_kms: f64 = parts[1].parse()?;
        let ev_kms: f64 = parts[2].parse()?;

        if r_kpc >= r_min && r_kpc <= r_max {
            rows.push(RcRow {
                r_kpc,
                vobs_kms,
                ev_kms,
            });
        }
    }

    if rows.is_empty() {
        bail!("no RC rows selected in range [{}, {}]", r_min, r_max);
    }

    Ok(rows)
}

fn parse_baryons(path: &str) -> Result<Vec<BaryonComponent>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    let mut rows = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line_no == 0 && line.starts_with("label,") {
            continue;
        }

        let parts: Vec<&str> = line.split(',').map(|v| v.trim()).collect();

        if parts.len() != 4 {
            bail!("bad baryon row {} in {}", line_no + 1, path);
        }

        rows.push(BaryonComponent {
            kind: parts[1].to_string(),
            mass_msun: parts[2].parse()?,
            scale_kpc: parts[3].parse()?,
        });
    }

    if rows.is_empty() {
        bail!("no baryon components loaded from {}", path);
    }

    Ok(rows)
}

fn baryon_component_velocity(component: &BaryonComponent, r_kpc: f64) -> Result<f64> {
    if r_kpc <= 0.0 {
        return Ok(0.0);
    }

    let menc = match component.kind.as_str() {
        "disk" | "exponential_disk" | "exponential-disk" => {
            let x = r_kpc / component.scale_kpc;
            component.mass_msun * (1.0 - (-x).exp() * (1.0 + x))
        }
        "bulge" | "hernquist" | "hernquist_bulge" | "hernquist-bulge" => {
            component.mass_msun * r_kpc.powi(2) / (r_kpc + component.scale_kpc).powi(2)
        }
        _ => bail!("unknown baryon kind: {}", component.kind),
    };

    Ok((G_KPC * menc / r_kpc).sqrt())
}

fn baryon_velocity_squared(baryons: &[BaryonComponent], r_kpc: f64) -> Result<f64> {
    let mut v2 = 0.0;

    for component in baryons {
        let v = baryon_component_velocity(component, r_kpc)?;
        v2 += v * v;
    }

    Ok(v2)
}

fn build_coeff_vectors(dim: usize, levels: i32) -> Vec<Vec<f64>> {
    let mut raw = Vec::new();
    let mut cur = vec![0i32; dim];

    fn rec(i: usize, dim: usize, levels: i32, cur: &mut [i32], raw: &mut Vec<Vec<i32>>) {
        if i == dim {
            if cur.iter().all(|&x| x == 0) {
                return;
            }

            // Remove global sign degeneracy: first nonzero coefficient is positive.
            for &x in cur.iter() {
                if x != 0 {
                    if x > 0 {
                        raw.push(cur.to_vec());
                    }
                    return;
                }
            }

            return;
        }

        for v in -levels..=levels {
            cur[i] = v;
            rec(i + 1, dim, levels, cur, raw);
        }
    }

    rec(0, dim, levels, &mut cur, &mut raw);

    let mut out = Vec::new();

    for v in raw {
        let norm = (v.iter().map(|&x| (x as f64) * (x as f64)).sum::<f64>()).sqrt();

        if norm <= 0.0 {
            continue;
        }

        out.push(v.iter().map(|&x| x as f64 / norm).collect());
    }

    out
}

fn coeff_tag(coeffs: &[f64]) -> String {
    coeffs
        .iter()
        .map(|c| format!("{:.6}", c))
        .collect::<Vec<_>>()
        .join(";")
}

fn score_candidate(
    rc: &[RcRow],
    baryons: &[BaryonComponent],
    state_ids: &[i32],
    coeffs: &[f64],
    a0: f64,
    dm_mass: f64,
    n_side: i32,
    extent_kpc: f64,
    softening_kpc: f64,
    rot: EulerRotation,
) -> Result<(f64, f64, f64)> {
    let terms: Vec<OrbitalTerm> = state_ids
        .iter()
        .zip(coeffs.iter())
        .map(|(&state_id, &coeff)| OrbitalTerm { state_id, coeff })
        .collect();

    let radii: Vec<f64> = rc.iter().map(|row| row.r_kpc).collect();
    let mut v_dm = vec![0.0f64; radii.len()];

    unsafe {
        compute_superposition_disk_curve_cuda(
            terms.as_ptr(),
            terms.len() as i32,
            a0,
            dm_mass,
            n_side,
            extent_kpc,
            softening_kpc,
            radii.as_ptr(),
            v_dm.as_mut_ptr(),
            radii.len(),
            rot,
        );
    }

    let mut sum_sq = 0.0;
    let mut chi2 = 0.0;

    for (i, row) in rc.iter().enumerate() {
        let vb2 = baryon_velocity_squared(baryons, row.r_kpc)?;
        let vtot = (v_dm[i] * v_dm[i] + vb2).sqrt();
        let residual = vtot - row.vobs_kms;

        sum_sq += residual * residual;

        if row.ev_kms > 0.0 {
            chi2 += (residual / row.ev_kms).powi(2);
        }
    }

    let n = rc.len() as f64;
    let rms = (sum_sq / n).sqrt();
    let chi2_per_point = chi2 / n;

    Ok((rms, chi2, chi2_per_point))
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 8 {
        bail!(
            "usage: scan_horb_cuda_mw_wavefunction <rc_csv> <baryons_csv> <out_csv> <r_min> <r_max> <n_side> <extent_kpc>"
        );
    }

    let rc_csv = &args[1];
    let baryons_csv = &args[2];
    let out_csv = &args[3];
    let r_min: f64 = args[4].parse()?;
    let r_max: f64 = args[5].parse()?;
    let n_side: i32 = args[6].parse()?;
    let extent_kpc: f64 = args[7].parse()?;

    let states = parse_list_string(
        &env::var("HORB_WAVE_STATE_LIST")
            .unwrap_or_else(|_| "3d_z2,3d_xy,3d_xz,3d_yz".to_string()),
    );

    let a0_values = parse_list_f64(
        &env::var("HORB_A0_LIST").unwrap_or_else(|_| "0.8,0.9,1.0,1.1,1.2,1.4".to_string()),
    )?;

    let mass_values = parse_list_f64(
        &env::var("HORB_DM_MASS_LIST").unwrap_or_else(|_| "1.5e11,2e11,2.5e11,3e11".to_string()),
    )?;

    let softening_values = parse_list_f64(
        &env::var("HORB_SOFTENING_LIST").unwrap_or_else(|_| "0.5".to_string()),
    )?;

    let coeff_levels: i32 = env::var("HORB_COEFF_LEVELS")
        .unwrap_or_else(|_| "2".to_string())
        .parse()?;

    let state_ids: Vec<i32> = states
        .iter()
        .map(|s| state_id(s))
        .collect::<Result<Vec<_>>>()?;

    let coeff_vectors = build_coeff_vectors(states.len(), coeff_levels);

    eprintln!("states={}", states.join(","));
    eprintln!("coeff_vectors={}", coeff_vectors.len());

    let rc = parse_standard_rc(rc_csv, r_min, r_max)?;
    let baryons = parse_baryons(baryons_csv)?;
    let rot = unsafe { make_identity_rotation_cuda() };

    let mut out = File::create(out_csv)?;

    writeln!(
        out,
        "states,coeffs,a0_star_kpc,dm_mass_msun,n_side,extent_kpc,softening_kpc,r_min_kpc,r_max_kpc,n,rms_kms,chi2,chi2_per_point"
    )?;

    for coeffs in coeff_vectors {
        for a0 in a0_values.iter().copied() {
            for dm_mass in mass_values.iter().copied() {
                for softening in softening_values.iter().copied() {
                    eprintln!(
                        "scoring psi=[{}] coeffs=[{}] a0={:.3} mass={:.6e} n={} eps={}",
                        states.join("+"),
                        coeff_tag(&coeffs),
                        a0,
                        dm_mass,
                        n_side,
                        softening
                    );

                    let (rms, chi2, chi2pt) = score_candidate(
                        &rc,
                        &baryons,
                        &state_ids,
                        &coeffs,
                        a0,
                        dm_mass,
                        n_side,
                        extent_kpc,
                        softening,
                        rot,
                    )?;

                    writeln!(
                        out,
                        "\"{}\",\"{}\",{:.6},{:.6e},{},{:.6},{:.6},{:.6},{:.6},{},{:.6},{:.6},{:.6}",
                        states.join(";"),
                        coeff_tag(&coeffs),
                        a0,
                        dm_mass,
                        n_side,
                        extent_kpc,
                        softening,
                        r_min,
                        r_max,
                        rc.len(),
                        rms,
                        chi2,
                        chi2pt
                    )?;

                    out.flush()?;
                }
            }
        }
    }

    Ok(())
}
