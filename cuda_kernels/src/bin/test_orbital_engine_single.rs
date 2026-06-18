use anyhow::{bail, Result};

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

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();

    let state = args.get(1).map(|s| s.as_str()).unwrap_or("3d_z2");
    let a0: f64 = args.get(2).map(|s| s.parse()).transpose()?.unwrap_or(0.9);
    let mass: f64 = args.get(3).map(|s| s.parse()).transpose()?.unwrap_or(2.0e11);
    let n_side: i32 = args.get(4).map(|s| s.parse()).transpose()?.unwrap_or(64);
    let extent: f64 = args.get(5).map(|s| s.parse()).transpose()?.unwrap_or(80.0);
    let softening: f64 = args.get(6).map(|s| s.parse()).transpose()?.unwrap_or(0.5);

    let radii: Vec<f64> = (1..=80).map(|i| i as f64).collect();
    let mut v_out = vec![0.0f64; radii.len()];

    let rot = unsafe { make_identity_rotation_cuda() };

    unsafe {
        compute_single_orbital_disk_curve_cuda(
            state_id(state)?,
            a0,
            mass,
            n_side,
            extent,
            softening,
            radii.as_ptr(),
            v_out.as_mut_ptr(),
            radii.len(),
            rot,
        );
    }

    println!("r_kpc,v_kms");
    for (r, v) in radii.iter().zip(v_out.iter()) {
        println!("{r:.6},{v:.6}");
    }

    Ok(())
}
