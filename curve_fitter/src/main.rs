// horb/curve_fitter/src/main.rs

use anyhow::{bail, Result};
use data_io::RotationCurve;
use orbital_basis::{
    BaryonComponent, BaryonComponentKind, BaryonicModel, ClassicalHalo, ClassicalHaloKind,
    DensityField, MultiBaryonicModel, OrbitalConfig,
};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        print_usage_and_exit()?;
    }

    let mode = &args[1];

    match mode.as_str() {
        "curve" => {
            if args.len() != 5 {
                print_usage_and_exit()?;
            }

            let field = build_field(&args[2], &args[3], &args[4])?;
            print_curve(&field)?;
        }

        "density" => {
            if args.len() != 5 {
                print_usage_and_exit()?;
            }

            let field = build_field(&args[2], &args[3], &args[4])?;
            print_density_axes(&field)?;
        }

        "xz" => {
            if args.len() != 5 {
                print_usage_and_exit()?;
            }

            let field = build_field(&args[2], &args[3], &args[4])?;
            print_xz_slice(&field)?;
        }

        "total" => {
            if args.len() != 9 {
                print_usage_and_exit()?;
            }

            let field = build_field(&args[2], &args[3], &args[4])?;

            let disk_mass: f64 = args[5].parse()?;
            let disk_scale: f64 = args[6].parse()?;
            let bulge_mass: f64 = args[7].parse()?;
            let bulge_scale: f64 = args[8].parse()?;

            let baryons = BaryonicModel::new(disk_mass, disk_scale, bulge_mass, bulge_scale)?;
            print_total_curve(&field, &baryons)?;
        }

        "halo" => {
            if args.len() != 6 {
                print_usage_and_exit()?;
            }

            let kind = ClassicalHaloKind::parse(&args[2])?;
            let scale_kpc: f64 = args[3].parse()?;
            let m_ref_msun: f64 = args[4].parse()?;
            let r_ref_kpc: f64 = args[5].parse()?;

            let halo = ClassicalHalo::from_mass_at_radius(kind, scale_kpc, m_ref_msun, r_ref_kpc)?;
            print_halo_curve(&halo)?;
        }

        "compare_dm" => {
            if args.len() != 9 {
                print_usage_and_exit()?;
            }

            let state = &args[2];
            let a0_star: f64 = args[3].parse()?;
            let dm_mass: f64 = args[4].parse()?;

            let scale_piso: f64 = args[5].parse()?;
            let scale_nfw: f64 = args[6].parse()?;
            let scale_burkert: f64 = args[7].parse()?;
            let r_ref_kpc: f64 = args[8].parse()?;

            let cfg = orbital_config(state, a0_star, dm_mass)?;
            let field = DensityField::new(cfg)?;

            let piso = ClassicalHalo::from_mass_at_radius(
                ClassicalHaloKind::PseudoIsothermal,
                scale_piso,
                dm_mass,
                r_ref_kpc,
            )?;

            let nfw = ClassicalHalo::from_mass_at_radius(
                ClassicalHaloKind::Nfw,
                scale_nfw,
                dm_mass,
                r_ref_kpc,
            )?;

            let burkert = ClassicalHalo::from_mass_at_radius(
                ClassicalHaloKind::Burkert,
                scale_burkert,
                dm_mass,
                r_ref_kpc,
            )?;

            print_dm_comparison_curve(&field, &piso, &nfw, &burkert)?;
        }

        "compare_dm_fixed_baselines" => {
            if args.len() != 10 {
                print_usage_and_exit()?;
            }

            let state = &args[2];
            let a0_star: f64 = args[3].parse()?;
            let horb_mass: f64 = args[4].parse()?;
            let baseline_mass: f64 = args[5].parse()?;

            let scale_piso: f64 = args[6].parse()?;
            let scale_nfw: f64 = args[7].parse()?;
            let scale_burkert: f64 = args[8].parse()?;
            let r_ref_kpc: f64 = args[9].parse()?;

            let cfg = orbital_config(state, a0_star, horb_mass)?;
            let field = DensityField::new(cfg)?;

            let piso = ClassicalHalo::from_mass_at_radius(
                ClassicalHaloKind::PseudoIsothermal,
                scale_piso,
                baseline_mass,
                r_ref_kpc,
            )?;

            let nfw = ClassicalHalo::from_mass_at_radius(
                ClassicalHaloKind::Nfw,
                scale_nfw,
                baseline_mass,
                r_ref_kpc,
            )?;

            let burkert = ClassicalHalo::from_mass_at_radius(
                ClassicalHaloKind::Burkert,
                scale_burkert,
                baseline_mass,
                r_ref_kpc,
            )?;

            print_dm_comparison_curve(&field, &piso, &nfw, &burkert)?;
        }
        "fit_standard_rc" => print_standard_rc_fit_score(&args[2..])?,

        "fit_standard_rc_multi_baryons" => print_standard_rc_fit_score_multi_baryons(&args[2..])?,
        _ => print_usage_and_exit()?,
    }

    Ok(())
}

fn print_usage_and_exit() -> Result<()> {
    bail!(
        "usage:\n\
         cargo run -p curve_fitter -- curve      <state> <a0_star_kpc> <dm_mass_msun>\n\
         cargo run -p curve_fitter -- density    <state> <a0_star_kpc> <dm_mass_msun>\n\
         cargo run -p curve_fitter -- xz         <state> <a0_star_kpc> <dm_mass_msun>\n\
         cargo run -p curve_fitter -- total      <state> <a0_star_kpc> <dm_mass_msun> <disk_mass_msun> <disk_scale_kpc> <bulge_mass_msun> <bulge_scale_kpc>\n\
         cargo run -p curve_fitter -- halo       <halo_kind> <scale_kpc> <m_ref_msun> <r_ref_kpc>\n\
         cargo run -p curve_fitter -- compare_dm <state> <a0_star_kpc> <dm_mass_msun> <piso_rc_kpc> <nfw_rs_kpc> <burkert_r0_kpc> <r_ref_kpc>\n\
         cargo run -p curve_fitter -- compare_dm_fixed_baselines <state> <a0_star_kpc> <horb_mass_msun> <baseline_mass_msun> <piso_rc_kpc> <nfw_rs_kpc> <burkert_r0_kpc> <r_ref_kpc>\n\
         \n\
         states: 1s, 3d_z2\n\
         halo_kind: piso, nfw, burkert\n\
         \n\
         examples:\n\
         cargo run -p curve_fitter -- halo nfw 15.0 1e11 80.0\n\
         cargo run -p curve_fitter -- compare_dm 3d_z2 1.5 1e11 5.0 15.0 8.0 80.0\n\
         cargo run -p curve_fitter -- compare_dm_fixed_baselines 3d_z2 1.5 5e10 1e11 5.0 15.0 8.0 80.0"
    )
}

fn build_field(state: &str, a0_star: &str, dm_mass: &str) -> Result<DensityField> {
    let a0_star: f64 = a0_star.parse()?;
    let dm_mass: f64 = dm_mass.parse()?;
    let cfg = orbital_config(state, a0_star, dm_mass)?;
    Ok(DensityField::new(cfg)?)
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

fn print_halo_curve(halo: &ClassicalHalo) -> Result<()> {
    println!("r_kpc,rho_Msun_per_kpc3,M_enc_Msun,v_circ_kms");

    let r_min = 0.05;
    let r_max = 80.0;
    let bins = 500;

    for i in 0..bins {
        let t = i as f64 / (bins - 1) as f64;
        let r = r_min + t * (r_max - r_min);

        let rho = halo.density(r);
        let m = halo.enclosed_mass(r);
        let v = halo.circular_velocity(r);

        println!("{:.6},{:.6e},{:.6e},{:.6}", r, rho, m, v);
    }

    Ok(())
}

fn print_dm_comparison_curve(
    horb: &DensityField,
    piso: &ClassicalHalo,
    nfw: &ClassicalHalo,
    burkert: &ClassicalHalo,
) -> Result<()> {
    println!("r_kpc,v_horb_kms,v_piso_kms,v_nfw_kms,v_burkert_kms,m_horb_msun,m_piso_msun,m_nfw_msun,m_burkert_msun");

    let r_min = 0.05;
    let r_max = 80.0;
    let bins = 500;
    let steps = 8000;

    for i in 0..bins {
        let t = i as f64 / (bins - 1) as f64;
        let r = r_min + t * (r_max - r_min);

        let m_horb = horb.enclosed_mass(r, steps)?;
        let v_horb = horb.circular_velocity_spherical(r, steps)?;

        let m_piso = piso.enclosed_mass(r);
        let v_piso = piso.circular_velocity(r);

        let m_nfw = nfw.enclosed_mass(r);
        let v_nfw = nfw.circular_velocity(r);

        let m_burkert = burkert.enclosed_mass(r);
        let v_burkert = burkert.circular_velocity(r);

        println!(
            "{:.6},{:.6},{:.6},{:.6},{:.6},{:.6e},{:.6e},{:.6e},{:.6e}",
            r, v_horb, v_piso, v_nfw, v_burkert, m_horb, m_piso, m_nfw, m_burkert
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

    let extent = std::env::var("HORB_XZ_EXTENT_KPC")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
        .unwrap_or(25.0);

    let n = std::env::var("HORB_XZ_N")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(201);

    if extent <= 0.0 {
        bail!("HORB_XZ_EXTENT_KPC must be positive");
    }

    if n < 3 {
        bail!("HORB_XZ_N must be >= 3");
    }

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

fn print_standard_rc_fit_score(args: &[String]) -> Result<()> {
    if args.len() != 10 {
        bail!(
            "usage: fit_standard_rc <csv> <state> <a0_star_kpc> <dm_mass_msun> \
             <disk_mass_msun> <disk_scale_kpc> <bulge_mass_msun> <bulge_scale_kpc> \
             <r_min_kpc> <r_max_kpc>"
        );
    }

    let csv_path = &args[0];
    let field = build_field(&args[1], &args[2], &args[3])?;

    let disk_mass = parse_f64(&args[4], "disk_mass_msun")?;
    let disk_scale = parse_f64(&args[5], "disk_scale_kpc")?;
    let bulge_mass = parse_f64(&args[6], "bulge_mass_msun")?;
    let bulge_scale = parse_f64(&args[7], "bulge_scale_kpc")?;
    let r_min = parse_f64(&args[8], "r_min_kpc")?;
    let r_max = parse_f64(&args[9], "r_max_kpc")?;

    let rc = RotationCurve::from_csv(csv_path)?;
    let baryons = BaryonicModel::new(disk_mass, disk_scale, bulge_mass, bulge_scale)?;

    let mut n = 0usize;
    let mut sum_sq = 0.0;
    let mut chi2 = 0.0;

    for row in rc.rows.iter() {
        if row.r_kpc < r_min || row.r_kpc > r_max {
            continue;
        }

        let v_dm = field.circular_velocity_spherical(row.r_kpc, 10_000)?;
        let v_total = if rc.has_baryons() {
            let vgas = row.vgas_kms.unwrap_or(0.0);
            let vdisk = row.vdisk_kms.unwrap_or(0.0);
            let vbul = row.vbul_kms.unwrap_or(0.0);

            (v_dm * v_dm + vgas * vgas + vdisk * vdisk + vbul * vbul).sqrt()
        } else {
            let v_disk = baryons.disk_velocity(row.r_kpc);
            let v_bulge = baryons.bulge_velocity(row.r_kpc);
            (v_dm * v_dm + v_disk * v_disk + v_bulge * v_bulge).sqrt()
        };

        let residual = v_total - row.vobs_kms;
        sum_sq += residual * residual;

        if row.ev_kms > 0.0 {
            chi2 += (residual / row.ev_kms).powi(2);
        }

        n += 1;
    }

    if n == 0 {
        bail!(
            "no rows selected from {} in range [{}, {}] kpc",
            csv_path,
            r_min,
            r_max
        );
    }

    let rms = (sum_sq / n as f64).sqrt();

    println!(
        "csv,state,a0_star_kpc,dm_mass_msun,disk_mass_msun,disk_scale_kpc,bulge_mass_msun,bulge_scale_kpc,r_min_kpc,r_max_kpc,n,rms_kms,chi2,chi2_per_point,has_baryons"
    );

    println!(
        "{},{},{:.6},{:.6e},{:.6e},{:.6},{:.6e},{:.6},{:.6},{:.6},{},{:.6},{:.6},{:.6},{}",
        csv_path,
        args[1],
        parse_f64(&args[2], "a0_star_kpc")?,
        parse_f64(&args[3], "dm_mass_msun")?,
        disk_mass,
        disk_scale,
        bulge_mass,
        bulge_scale,
        r_min,
        r_max,
        n,
        rms,
        chi2,
        chi2 / n as f64,
        rc.has_baryons()
    );

    Ok(())
}

fn parse_f64(value: &str, name: &str) -> Result<f64> {
    value
        .parse::<f64>()
        .map_err(|err| anyhow::anyhow!("failed to parse {}='{}' as f64: {}", name, value, err))
}

fn load_multi_baryons_csv(path: &str) -> Result<MultiBaryonicModel> {
    let content = std::fs::read_to_string(path)?;

    let mut components = Vec::new();

    for (line_no, line) in content.lines().enumerate() {
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line_no == 0 && line.to_lowercase().starts_with("label,") {
            continue;
        }

        let parts: Vec<&str> = line.split(',').map(|v| v.trim()).collect();

        if parts.len() != 4 {
            anyhow::bail!(
                "bad baryon CSV row {} in {}: expected 4 columns, got {}",
                line_no + 1,
                path,
                parts.len()
            );
        }

        let label = parts[0].to_string();

        let kind = BaryonComponentKind::parse(parts[1]).ok_or_else(|| {
            anyhow::anyhow!(
                "bad baryon component kind '{}' on row {} in {}",
                parts[1],
                line_no + 1,
                path
            )
        })?;

        let mass_msun = parse_f64(parts[2], "mass_msun")?;
        let scale_kpc = parse_f64(parts[3], "scale_kpc")?;

        components.push(BaryonComponent::new(label, kind, mass_msun, scale_kpc)?);
    }

    MultiBaryonicModel::new(components)
}

fn print_standard_rc_fit_score_multi_baryons(args: &[String]) -> Result<()> {
    if args.len() != 7 {
        bail!(
            "usage: fit_standard_rc_multi_baryons <rc_csv> <baryons_csv> <state> \
             <a0_star_kpc> <dm_mass_msun> <r_min_kpc> <r_max_kpc>"
        );
    }

    let rc_csv = &args[0];
    let baryons_csv = &args[1];
    let state = &args[2];
    let a0_star = &args[3];
    let dm_mass = &args[4];
    let r_min = parse_f64(&args[5], "r_min_kpc")?;
    let r_max = parse_f64(&args[6], "r_max_kpc")?;

    let field = build_field(state, a0_star, dm_mass)?;
    let rc = RotationCurve::from_csv(rc_csv)?;
    let baryons = load_multi_baryons_csv(baryons_csv)?;

    let mut n = 0usize;
    let mut sum_sq = 0.0;
    let mut chi2 = 0.0;

    for row in rc.rows.iter() {
        if row.r_kpc < r_min || row.r_kpc > r_max {
            continue;
        }

        let v_dm = field.circular_velocity_spherical(row.r_kpc, 10_000)?;

        let v_total = if rc.has_baryons() {
            let vgas = row.vgas_kms.unwrap_or(0.0);
            let vdisk = row.vdisk_kms.unwrap_or(0.0);
            let vbul = row.vbul_kms.unwrap_or(0.0);

            (v_dm * v_dm + vgas * vgas + vdisk * vdisk + vbul * vbul).sqrt()
        } else {
            baryons.total_velocity(row.r_kpc, v_dm)
        };

        let residual = v_total - row.vobs_kms;

        sum_sq += residual * residual;

        if row.ev_kms > 0.0 {
            chi2 += (residual / row.ev_kms).powi(2);
        }

        n += 1;
    }

    if n == 0 {
        bail!(
            "no rows selected from {} in range [{}, {}] kpc",
            rc_csv,
            r_min,
            r_max
        );
    }

    let rms = (sum_sq / n as f64).sqrt();

    println!(
        "rc_csv,baryons_csv,state,a0_star_kpc,dm_mass_msun,r_min_kpc,r_max_kpc,n,rms_kms,chi2,chi2_per_point,has_rc_baryons,n_baryon_components"
    );

    println!(
        "{},{},{},{:.6},{:.6e},{:.6},{:.6},{},{:.6},{:.6},{:.6},{},{}",
        rc_csv,
        baryons_csv,
        state,
        parse_f64(a0_star, "a0_star_kpc")?,
        parse_f64(dm_mass, "dm_mass_msun")?,
        r_min,
        r_max,
        n,
        rms,
        chi2,
        chi2 / n as f64,
        rc.has_baryons(),
        baryons.components.len()
    );

    Ok(())
}
