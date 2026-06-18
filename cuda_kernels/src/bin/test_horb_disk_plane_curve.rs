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

fn rho_3d_z2_unnormalized(x: f64, y: f64, z: f64, a0: f64) -> f64 {
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();

    if r <= 1.0e-12 {
        return 0.0;
    }

    let cos_theta = z / r;
    let angular = 3.0 * cos_theta * cos_theta - 1.0;

    // Same morphology density used in the CPU field:
    // rho ∝ r^4 exp(-2r / 3a0) (3cos²θ - 1)²
    let radial = r2 * r2 * (-2.0 * r / (3.0 * a0)).exp();

    radial * angular * angular
}

fn build_mass_cells(
    n_side: usize,
    extent_kpc: f64,
    a0: f64,
    total_mass_msun: f64,
) -> Result<Vec<MassCell>> {
    if n_side < 8 {
        bail!("n_side must be >= 8");
    }

    if extent_kpc <= 0.0 {
        bail!("extent_kpc must be positive");
    }

    if a0 <= 0.0 {
        bail!("a0 must be positive");
    }

    if total_mass_msun <= 0.0 {
        bail!("total_mass_msun must be positive");
    }

    let dx = 2.0 * extent_kpc / n_side as f64;

    let mut raw = Vec::with_capacity(n_side * n_side * n_side);
    let mut rho_sum = 0.0f64;

    for iz in 0..n_side {
        let z = -extent_kpc + (iz as f64 + 0.5) * dx;

        for iy in 0..n_side {
            let y = -extent_kpc + (iy as f64 + 0.5) * dx;

            for ix in 0..n_side {
                let x = -extent_kpc + (ix as f64 + 0.5) * dx;

                let rho = rho_3d_z2_unnormalized(x, y, z, a0);
                rho_sum += rho;

                raw.push((x, y, z, rho));
            }
        }
    }

    if rho_sum <= 0.0 {
        bail!("rho_sum was non-positive");
    }

    let cells = raw
        .into_iter()
        .map(|(x, y, z, rho)| MassCell {
            x,
            y,
            z,
            m: total_mass_msun * rho / rho_sum,
        })
        .collect();

    Ok(cells)
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    let a0: f64 = args.get(1).map(|s| s.parse()).transpose()?.unwrap_or(1.0);
    let total_mass: f64 = args
        .get(2)
        .map(|s| s.parse())
        .transpose()?
        .unwrap_or(2.5e11);
    let n_side: usize = args.get(3).map(|s| s.parse()).transpose()?.unwrap_or(96);
    let extent_kpc: f64 = args.get(4).map(|s| s.parse()).transpose()?.unwrap_or(80.0);
    let softening_kpc: f64 = args.get(5).map(|s| s.parse()).transpose()?.unwrap_or(0.25);

    eprintln!("building mass grid...");
    eprintln!("a0_kpc={a0}");
    eprintln!("total_mass_msun={total_mass:.6e}");
    eprintln!("n_side={n_side}");
    eprintln!("extent_kpc={extent_kpc}");
    eprintln!("softening_kpc={softening_kpc}");

    let cells = build_mass_cells(n_side, extent_kpc, a0, total_mass)?;
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

    let out_path = format!(
        "curves/cuda_horb_diskplane_a{:.3}_m{:.3e}_n{}.csv",
        a0, total_mass, n_side
    );

    std::fs::create_dir_all("curves")?;

    let mut f = File::create(&out_path)?;
    writeln!(f, "r_kpc,v_horb_cuda_kms")?;

    for (r, v) in radii.iter().zip(v_out.iter()) {
        writeln!(f, "{r:.6},{v:.6}")?;
    }

    eprintln!("wrote {out_path}");

    println!("r_kpc,v_horb_cuda_kms");
    for (r, v) in radii.iter().zip(v_out.iter()).take(20) {
        println!("{r:.6},{v:.6}");
    }

    Ok(())
}
