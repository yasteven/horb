use anyhow::{bail, Result};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct MassCell {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub m: f64,
}

#[link(name = "horb_cuda", kind = "dylib")]
unsafe extern "C" {
    fn compute_disk_plane_rotation_curve(
        cells: *const MassCell,
        cell_count: usize,
        radii: *const f64,
        v_out: *mut f64,
        radius_count: usize,
        softening_kpc: f64,
    );
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

fn parse_list_f64(value: &str) -> Result<Vec<f64>> {
    value
        .split(',')
        .map(|v| Ok(v.trim().parse::<f64>()?))
        .collect()
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

fn rho_hydrogenic_unnormalized(state: &str, x: f64, y: f64, z: f64, a0: f64) -> Result<f64> {
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();

    if r <= 1.0e-12 {
        return Ok(0.0);
    }

    let rho = match state {
        // 1s: |R_10 Y_00|^2, ignoring constants.
        "1s" => (-2.0 * r / a0).exp(),

        // 2p real states: density ∝ r^2 exp(-r/a0) angular^2.
        "2p_z" => {
            let angular = z / r;
            r2 * (-r / a0).exp() * angular * angular
        }
        "2p_x" => {
            let angular = x / r;
            r2 * (-r / a0).exp() * angular * angular
        }
        "2p_y" => {
            let angular = y / r;
            r2 * (-r / a0).exp() * angular * angular
        }

        // 3d real states: density ∝ r^4 exp(-2r/3a0) angular^2.
        "3d_z2" => {
            let cos_theta = z / r;
            let angular = 3.0 * cos_theta * cos_theta - 1.0;
            r2 * r2 * (-2.0 * r / (3.0 * a0)).exp() * angular * angular
        }
        "3d_x2_y2" => {
            let angular = (x * x - y * y) / r2;
            r2 * r2 * (-2.0 * r / (3.0 * a0)).exp() * angular * angular
        }
        "3d_xy" => {
            let angular = x * y / r2;
            r2 * r2 * (-2.0 * r / (3.0 * a0)).exp() * angular * angular
        }
        "3d_xz" => {
            let angular = x * z / r2;
            r2 * r2 * (-2.0 * r / (3.0 * a0)).exp() * angular * angular
        }
        "3d_yz" => {
            let angular = y * z / r2;
            r2 * r2 * (-2.0 * r / (3.0 * a0)).exp() * angular * angular
        }

        _ => bail!("unknown orbital state: {}", state),
    };

    Ok(rho)
}

fn build_mass_cells(
    state: &str,
    n_side: usize,
    extent_kpc: f64,
    a0: f64,
    total_mass_msun: f64,
) -> Result<Vec<MassCell>> {
    if n_side < 8 {
        bail!("n_side must be >= 8");
    }

    let dx = 2.0 * extent_kpc / n_side as f64;

    let mut raw = Vec::with_capacity(n_side * n_side * n_side);
    let mut rho_sum = 0.0;

    for iz in 0..n_side {
        let z = -extent_kpc + (iz as f64 + 0.5) * dx;

        for iy in 0..n_side {
            let y = -extent_kpc + (iy as f64 + 0.5) * dx;

            for ix in 0..n_side {
                let x = -extent_kpc + (ix as f64 + 0.5) * dx;
                let rho = rho_hydrogenic_unnormalized(state, x, y, z, a0)?;

                rho_sum += rho;
                raw.push((x, y, z, rho));
            }
        }
    }

    if rho_sum <= 0.0 {
        bail!("rho_sum was non-positive for state {}", state);
    }

    Ok(raw
        .into_iter()
        .map(|(x, y, z, rho)| MassCell {
            x,
            y,
            z,
            m: total_mass_msun * rho / rho_sum,
        })
        .collect())
}

fn score_candidate(
    rc: &[RcRow],
    baryons: &[BaryonComponent],
    state: &str,
    a0: f64,
    dm_mass: f64,
    n_side: usize,
    extent_kpc: f64,
    softening_kpc: f64,
) -> Result<(f64, f64, f64)> {
    let cells = build_mass_cells(state, n_side, extent_kpc, a0, dm_mass)?;

    let radii: Vec<f64> = rc.iter().map(|row| row.r_kpc).collect();
    let mut v_dm = vec![0.0f64; radii.len()];

    unsafe {
        compute_disk_plane_rotation_curve(
            cells.as_ptr(),
            cells.len(),
            radii.as_ptr(),
            v_dm.as_mut_ptr(),
            radii.len(),
            softening_kpc,
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
            "usage: scan_horb_cuda_mw_basis <rc_csv> <baryons_csv> <out_csv> <r_min> <r_max> <n_side> <extent_kpc>"
        );
    }

    let rc_csv = &args[1];
    let baryons_csv = &args[2];
    let out_csv = &args[3];
    let r_min: f64 = args[4].parse()?;
    let r_max: f64 = args[5].parse()?;
    let n_side: usize = args[6].parse()?;
    let extent_kpc: f64 = args[7].parse()?;

    let state_list = env::var("HORB_STATE_LIST")
        .unwrap_or_else(|_| "1s,2p_z,2p_x,2p_y,3d_z2,3d_x2_y2,3d_xy,3d_xz,3d_yz".to_string());
    let a0_list = env::var("HORB_A0_LIST")
        .unwrap_or_else(|_| "0.6,0.8,0.9,1.0,1.1,1.2,1.4".to_string());
    let mass_list = env::var("HORB_DM_MASS_LIST")
        .unwrap_or_else(|_| "1.5e11,2e11,2.5e11,3e11".to_string());
    let softening_list = env::var("HORB_SOFTENING_LIST")
        .unwrap_or_else(|_| "0.25,0.5".to_string());

    let states = parse_list_string(&state_list);
    let a0_values = parse_list_f64(&a0_list)?;
    let mass_values = parse_list_f64(&mass_list)?;
    let softening_values = parse_list_f64(&softening_list)?;

    let rc = parse_standard_rc(rc_csv, r_min, r_max)?;
    let baryons = parse_baryons(baryons_csv)?;

    let mut out = File::create(out_csv)?;

    writeln!(
        out,
        "state,a0_star_kpc,dm_mass_msun,n_side,extent_kpc,softening_kpc,r_min_kpc,r_max_kpc,n,rms_kms,chi2,chi2_per_point"
    )?;

    for state in states {
        for a0 in a0_values.iter().copied() {
            for dm_mass in mass_values.iter().copied() {
                for softening in softening_values.iter().copied() {
                    eprintln!(
                        "scoring {} a0={:.3} mass={:.6e} n_side={} extent={} softening={}",
                        state, a0, dm_mass, n_side, extent_kpc, softening
                    );

                    let (rms, chi2, chi2pt) = score_candidate(
                        &rc,
                        &baryons,
                        &state,
                        a0,
                        dm_mass,
                        n_side,
                        extent_kpc,
                        softening,
                    )?;

                    writeln!(
                        out,
                        "{},{:.6},{:.6e},{},{:.6},{:.6},{:.6},{:.6},{},{:.6},{:.6},{:.6}",
                        state,
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
