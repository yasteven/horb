use anyhow::{bail, Result};
use std::env;
use std::fs::File;
use std::io::Write;

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

fn rho_hydrogenic_unnormalized(state: &str, x: f64, y: f64, z: f64, a0: f64) -> Result<f64> {
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();

    if r <= 1.0e-12 {
        return Ok(0.0);
    }

    let rho = match state {
        "1s" => (-2.0 * r / a0).exp(),

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

fn safe_state_tag(state: &str) -> String {
    state.replace('/', "_")
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 7 {
        bail!(
            "usage: write_horb_cuda_curve_basis <state> <a0_kpc> <mass_msun> <n_side> <extent_kpc> <softening_kpc>"
        );
    }

    let state = &args[1];
    let a0: f64 = args[2].parse()?;
    let total_mass: f64 = args[3].parse()?;
    let n_side: usize = args[4].parse()?;
    let extent_kpc: f64 = args[5].parse()?;
    let softening_kpc: f64 = args[6].parse()?;

    eprintln!("building mass grid...");
    eprintln!("state={state}");
    eprintln!("a0_kpc={a0}");
    eprintln!("total_mass_msun={total_mass:.6e}");
    eprintln!("n_side={n_side}");
    eprintln!("extent_kpc={extent_kpc}");
    eprintln!("softening_kpc={softening_kpc}");

    let cells = build_mass_cells(state, n_side, extent_kpc, a0, total_mass)?;
    let mass_sum: f64 = cells.iter().map(|c| c.m).sum();

    eprintln!("cells={}", cells.len());
    eprintln!("mass_sum={mass_sum:.6e}");

    let radii: Vec<f64> = (1..=160).map(|i| i as f64 * 0.5).collect();
    let mut v_out = vec![0.0f64; radii.len()];

    unsafe {
        compute_disk_plane_rotation_curve(
            cells.as_ptr(),
            cells.len(),
            radii.as_ptr(),
            v_out.as_mut_ptr(),
            radii.len(),
            softening_kpc,
        );
    }

    std::fs::create_dir_all("curves")?;

    let out_path = format!(
        "curves/cuda_horb_basis_{}_a{:.3}_m{:.3e}_n{}_eps{:.3}.csv",
        safe_state_tag(state),
        a0,
        total_mass,
        n_side,
        softening_kpc
    )
    .replace("e+", "e");

    let mut f = File::create(&out_path)?;
    writeln!(f, "r_kpc,v_horb_cuda_kms")?;

    for (r, v) in radii.iter().zip(v_out.iter()) {
        writeln!(f, "{r:.6},{v:.6}")?;
    }

    eprintln!("wrote {out_path}");
    println!("{out_path}");

    Ok(())
}
