// horb/curve_fitter/src/main.rs

use anyhow::{bail, Result};
use orbital_basis::{BaryonicModel, DensityField, OrbitalConfig};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 5 {
        print_usage_and_exit()?;
    }

    let mode = &args[1];
    let state = &args[2];
    let a0_star: f64 = args[3].parse()?;
    let dm_mass: f64 = args[4].parse()?;

    let cfg = orbital_config(state, a0_star, dm_mass)?;
    let field = DensityField::new(cfg)?;

    match mode.as_str() {
        "curve" => {
            if args.len() != 5 {
                print_usage_and_exit()?;
            }
            print_curve(&field)?;
        }
        "density" => {
            if args.len() != 5 {
                print_usage_and_exit()?;
            }
            print_density_axes(&field)?;
        }
        "xz" => {
            if args.len() != 5 {
                print_usage_and_exit()?;
            }
            print_xz_slice(&field)?;
        }
        "total" => {
            if args.len() != 9 {
                print_usage_and_exit()?;
            }

            let disk_mass: f64 = args[5].parse()?;
            let disk_scale: f64 = args[6].parse()?;
            let bulge_mass: f64 = args[7].parse()?;
            let bulge_scale: f64 = args[8].parse()?;

            let baryons = BaryonicModel::new(disk_mass, disk_scale, bulge_mass, bulge_scale)?;
            print_total_curve(&field, &baryons)?;
        }
        _ => bail!("unknown mode '{}'", mode),
    }

    Ok(())
}

fn print_usage_and_exit() -> Result<()> {
    bail!(
        "usage:\n\
         cargo run -p curve_fitter -- curve   <state> <a0_star_kpc> <dm_mass_msun>\n\
         cargo run -p curve_fitter -- density <state> <a0_star_kpc> <dm_mass_msun>\n\
         cargo run -p curve_fitter -- xz      <state> <a0_star_kpc> <dm_mass_msun>\n\
         cargo run -p curve_fitter -- total   <state> <a0_star_kpc> <dm_mass_msun> <disk_mass_msun> <disk_scale_kpc> <bulge_mass_msun> <bulge_scale_kpc>\n\
         \n\
         states: 1s, 3d_z2\n\
         example:\n\
         cargo run -p curve_fitter -- total 3d_z2 1.5 1e11 6e10 3.0 1e10 0.7"
    )
}

fn orbital_config(state: &str, a0_star: f64, dm_mass: f64) -> Result<OrbitalConfig> {
    let cfg = match state {
        "1s" => OrbitalConfig::ground_state(a0_star, dm_mass)?,
        "3d_z2" => OrbitalConfig::d_z2(a0_star, dm_mass)?,
        _ => bail!("unknown state '{}'; supported: 1s, 3d_z2", state),
    };

    Ok(cfg)
}

fn print_curve(field: &DensityField) -> Result<()> {
    println!("r_kpc,M_enc_Msun,v_circ_kms");

    let r_min = 0.05;
    let r_max = 80.0;
    let bins = 500;
    let steps = 8000;

    for i in 0..bins {
        let t = i as f64 / (bins - 1) as f64;
        let r = r_min + t * (r_max - r_min);

        let m_enc = field.enclosed_mass(r, steps)?;
        let v = field.circular_velocity_spherical(r, steps)?;

        println!("{:.6},{:.6e},{:.6}", r, m_enc, v);
    }

    Ok(())
}

fn print_total_curve(field: &DensityField, baryons: &BaryonicModel) -> Result<()> {
    println!("r_kpc,M_dm_enc_Msun,v_dm_kms,v_disk_kms,v_bulge_kms,v_baryon_kms,v_total_kms");

    let r_min = 0.05;
    let r_max = 80.0;
    let bins = 500;
    let steps = 8000;

    for i in 0..bins {
        let t = i as f64 / (bins - 1) as f64;
        let r = r_min + t * (r_max - r_min);

        let m_dm = field.enclosed_mass(r, steps)?;
        let v_dm = field.circular_velocity_spherical(r, steps)?;
        let v_disk = baryons.disk_velocity(r);
        let v_bulge = baryons.bulge_velocity(r);
        let v_baryon = baryons.baryon_velocity(r);
        let v_total = BaryonicModel::total_velocity(v_dm, v_disk, v_bulge);

        println!(
            "{:.6},{:.6e},{:.6},{:.6},{:.6},{:.6},{:.6}",
            r, m_dm, v_dm, v_disk, v_bulge, v_baryon, v_total
        );
    }

    Ok(())
}

fn print_density_axes(field: &DensityField) -> Result<()> {
    println!("r_kpc,rho_pole_Msun_per_kpc3,rho_equator_Msun_per_kpc3,pole_to_equator_ratio");

    let r_min = 0.05;
    let r_max = 40.0;
    let bins = 500;

    for i in 0..bins {
        let t = i as f64 / (bins - 1) as f64;
        let r = r_min + t * (r_max - r_min);

        let rho_pole = field.rho(r, 0.0, 0.0)?;
        let rho_equator = field.rho(r, std::f64::consts::FRAC_PI_2, 0.0)?;
        let ratio = rho_pole / rho_equator;

        println!("{:.6},{:.6e},{:.6e},{:.6}", r, rho_pole, rho_equator, ratio);
    }

    Ok(())
}

fn print_xz_slice(field: &DensityField) -> Result<()> {
    println!("x_kpc,z_kpc,r_kpc,theta_rad,rho_Msun_per_kpc3");

    let extent = 25.0;
    let n = 201;

    for iz in 0..n {
        let z = -extent + 2.0 * extent * iz as f64 / (n - 1) as f64;

        for ix in 0..n {
            let x = -extent + 2.0 * extent * ix as f64 / (n - 1) as f64;

            let r = (x * x + z * z).sqrt();

            let theta = if r == 0.0 { 0.0 } else { (z / r).acos() };

            let phi = 0.0;
            let rho = field.rho(r, theta, phi)?;

            println!("{:.6},{:.6},{:.6},{:.6},{:.6e}", x, z, r, theta, rho);
        }
    }

    Ok(())
}
