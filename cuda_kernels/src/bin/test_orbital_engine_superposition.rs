use anyhow::Result;

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

fn main() -> Result<()> {
    // Example real superposition:
    // psi = 1.0 * 3d_z2 + 0.25 * 3d_yz
    //
    // state IDs:
    // 20 = 3d_z2
    // 24 = 3d_yz
    let terms = vec![
        OrbitalTerm {
            state_id: 20,
            coeff: 1.0,
        },
        OrbitalTerm {
            state_id: 24,
            coeff: 0.25,
        },
    ];

    let a0 = 0.9;
    let mass = 2.0e11;
    let n_side = 64;
    let extent = 80.0;
    let softening = 0.5;

    let radii: Vec<f64> = (1..=80).map(|i| i as f64).collect();
    let mut v_out = vec![0.0f64; radii.len()];

    let rot = unsafe { make_identity_rotation_cuda() };

    unsafe {
        compute_superposition_disk_curve_cuda(
            terms.as_ptr(),
            terms.len() as i32,
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
