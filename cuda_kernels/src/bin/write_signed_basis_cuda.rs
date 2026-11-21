use anyhow::{bail, Result};
use orbital_basis::signed_gravity::{mass_cells, MassCell};

#[link(name = "horb_cuda", kind = "dylib")]
unsafe extern "C" {
    fn compute_disk_plane_signed_v2_cuda(
        cells: *const MassCell,
        count: usize,
        radii: *const f64,
        out: *mut f64,
        radii_count: usize,
        softening: f64,
    ) -> i32;
}

fn main() -> Result<()> {
    data_io::signed_basis::write_library(
        &std::env::args().skip(1).collect::<Vec<_>>(),
        "cuda_quadrature",
        |cfg, radii, options| {
            let cells = mass_cells(cfg, options.nr, options.nmu, options.nphi)?;
            let mass: f64 = cells.iter().map(|c| c.m).sum();
            if (mass / cfg.dm_mass - 1.0).abs() > 1e-4 {
                bail!("quadrature mass differs from analytic normalization by more than 1e-4; increase --n-radial/--n-mu");
            }
            let mut out = vec![f64::NAN; radii.len()];
            let status = unsafe {
                compute_disk_plane_signed_v2_cuda(
                    cells.as_ptr(),
                    cells.len(),
                    radii.as_ptr(),
                    out.as_mut_ptr(),
                    radii.len(),
                    options.softening,
                )
            };
            if status != 1 {
                bail!("CUDA signed-force calculation failed");
            }
            Ok(out)
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cuda_preserves_signed_cells_and_outward_force() {
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
            MassCell {
                x: 0.0,
                y: 3.0,
                z: 2.0,
                m: -1e9,
            },
        ];
        let radii = [1.0, 15.0, 30.0];
        let mut out = [f64::NAN; 3];
        assert_eq!(
            unsafe {
                compute_disk_plane_signed_v2_cuda(
                    cells.as_ptr(),
                    cells.len(),
                    radii.as_ptr(),
                    out.as_mut_ptr(),
                    radii.len(),
                    0.1,
                )
            },
            1
        );
        for (r, k) in radii.iter().zip(out) {
            let cpu = orbital_basis::signed_gravity::direct_v2(&cells, *r, 0.1);
            assert!((k - cpu).abs() < 1e-9 * cpu.abs().max(1.0));
        }
        assert!(out[0] < 0.0);
    }
}
