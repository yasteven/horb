use anyhow::{bail, Result};
use std::env;
use std::fs::File;
use std::io::Write;

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
    fn compute_single_orbital_disk_curve_cuda(
        state_id: i32,
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

fn parse_list_f64(s: &str) -> Result<Vec<f64>> {
    s.split(',')
        .map(|x| Ok(x.trim().parse::<f64>()?))
        .collect()
}

fn parse_list_string(s: &str) -> Vec<String> {
    s.split(',')
        .map(|x| x.trim().to_string())
        .filter(|x| !x.is_empty())
        .collect()
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 6 {
        bail!(
            "usage: write_cuda_orbital_basis_library <out_csv> <n_side> <extent_kpc> <softening_kpc> <ref_mass_msun>"
        );
    }

    let out_csv = &args[1];
    let n_side: i32 = args[2].parse()?;
    let extent_kpc: f64 = args[3].parse()?;
    let softening_kpc: f64 = args[4].parse()?;
    let ref_mass_msun: f64 = args[5].parse()?;

    let states = parse_list_string(
        &env::var("HORB_STATE_LIST")
            .unwrap_or_else(|_| "1s,2p_z,2p_x,2p_y,3d_z2,3d_x2_y2,3d_xy,3d_xz,3d_yz".to_string()),
    );

    let a0_values = parse_list_f64(
        &env::var("HORB_A0_LIST")
            .unwrap_or_else(|_| "0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4".to_string()),
    )?;

    let radii: Vec<f64> = (1..=160).map(|i| i as f64 * 0.5).collect();
    let rot = unsafe { make_identity_rotation_cuda() };

    std::fs::create_dir_all("reports")?;
    std::fs::create_dir_all("curves")?;

    let mut out = File::create(out_csv)?;

    writeln!(
        out,
        "basis_id,state,a0_star_kpc,ref_mass_msun,n_side,extent_kpc,softening_kpc,r_kpc,v_cuda_kms"
    )?;

    let mut basis_id = 0usize;

    for state in states {
        for a0 in a0_values.iter().copied() {
            eprintln!(
                "basis_id={} state={} a0={} ref_mass={:.6e} n_side={} extent={} softening={}",
                basis_id, state, a0, ref_mass_msun, n_side, extent_kpc, softening_kpc
            );

            let mut v_out = vec![0.0f64; radii.len()];

            unsafe {
                compute_single_orbital_disk_curve_cuda(
                    state_id(&state)?,
                    a0,
                    ref_mass_msun,
                    n_side,
                    extent_kpc,
                    softening_kpc,
                    radii.as_ptr(),
                    v_out.as_mut_ptr(),
                    radii.len(),
                    rot,
                );
            }

            for (r, v) in radii.iter().zip(v_out.iter()) {
                writeln!(
                    out,
                    "{},{},{:.6},{:.6e},{},{:.6},{:.6},{:.6},{:.6}",
                    basis_id,
                    state,
                    a0,
                    ref_mass_msun,
                    n_side,
                    extent_kpc,
                    softening_kpc,
                    r,
                    v
                )?;
            }

            out.flush()?;
            basis_id += 1;
        }
    }

    eprintln!("wrote {}", out_csv);

    Ok(())
}
