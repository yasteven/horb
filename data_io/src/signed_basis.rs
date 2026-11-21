//! Versioned, signed force-basis export shared by CPU and CUDA entry points.
use anyhow::{bail, Context, Result};
use orbital_basis::{signed_gravity::state_config, OrbitalConfig, Superposition};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
};

pub const DEFAULT_STATES: &str = "1s,2s,2p_z,3s,3p_z,3d_z2";
pub const DEFAULT_SCALES: &str = "0.5,1,2,4,8";

pub fn options(args: &[String], allowed: &[&str]) -> Result<BTreeMap<String, String>> {
    if args.len() % 2 != 0 {
        bail!("options must be --name VALUE pairs");
    }
    let mut values = BTreeMap::new();
    for pair in args.chunks(2) {
        if !allowed.contains(&pair[0].as_str())
            || values.insert(pair[0].clone(), pair[1].clone()).is_some()
        {
            bail!("unknown or duplicate option: {}", pair[0]);
        }
    }
    Ok(values)
}

fn value<'a>(o: &'a BTreeMap<String, String>, name: &str, default: &'a str) -> &'a str {
    o.get(name).map(String::as_str).unwrap_or(default)
}

fn number(o: &BTreeMap<String, String>, name: &str, default: &str) -> Result<f64> {
    let x: f64 = value(o, name, default)
        .parse()
        .with_context(|| format!("invalid {name}"))?;
    if !x.is_finite() {
        bail!("{name} must be finite");
    }
    Ok(x)
}

pub struct BasisOptions {
    pub nr: usize,
    pub nmu: usize,
    pub nphi: usize,
    pub softening: f64,
}

#[derive(Serialize)]
struct BasisRow<'a> {
    schema_version: u32,
    basis_id: usize,
    state: &'a str,
    n: u32,
    l: u32,
    m: i32,
    a0_star_kpc: f64,
    ref_mass_msun: f64,
    gravity_backend: &'a str,
    geometry: &'a str,
    normalization: &'a str,
    n_radial: usize,
    n_mu: usize,
    n_phi: usize,
    softening_kpc: f64,
    r_kpc: f64,
    v2_signed_kms2: f64,
}

fn output_file(path: &Path) -> Result<tempfile::NamedTempFile> {
    let parent = path
        .parent()
        .filter(|x| !x.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent)?;
    Ok(tempfile::NamedTempFile::new_in(parent)?)
}

pub fn write_library<F>(args: &[String], backend: &str, mut calculate: F) -> Result<()>
where
    F: FnMut(&OrbitalConfig, &[f64], &BasisOptions) -> Result<Vec<f64>>,
{
    if args == ["--help"] {
        println!("signed_basis --out FILE [--states {DEFAULT_STATES}] [--scales {DEFAULT_SCALES}] [--ref-mass 1e11] [--r-min 0.5 --r-max 100 --r-step 0.5] [--rc-csv FILE] [--n-radial 128 --n-mu 32 --n-phi 64 --softening 0]");
        return Ok(());
    }
    let o = options(
        args,
        &[
            "--out",
            "--states",
            "--scales",
            "--ref-mass",
            "--r-min",
            "--r-max",
            "--r-step",
            "--rc-csv",
            "--n-radial",
            "--n-mu",
            "--n-phi",
            "--softening",
        ],
    )?;
    let out = PathBuf::from(o.get("--out").context("--out is required")?);
    let opts = BasisOptions {
        nr: value(&o, "--n-radial", "128").parse()?,
        nmu: value(&o, "--n-mu", "32").parse()?,
        nphi: value(&o, "--n-phi", "64").parse()?,
        softening: number(&o, "--softening", "0")?,
    };
    for order in [opts.nr, opts.nmu, opts.nphi] {
        if !(8..=512).contains(&order) {
            bail!("quadrature orders must be 8..512");
        }
    }
    if opts.softening < 0.0 {
        bail!("softening must be nonnegative");
    }
    if backend == "cpu_multipole" && opts.softening != 0.0 {
        bail!("CPU multipoles require --softening 0");
    }
    let lo = number(&o, "--r-min", "0.5")?;
    let hi = number(&o, "--r-max", "100")?;
    let step = number(&o, "--r-step", "0.5")?;
    if !(0.0 < lo && lo < hi && step > 0.0) || (hi - lo) / step > 100_000.0 {
        bail!("invalid radius grid (at most 100000 intervals)");
    }
    let mut radii: Vec<f64> = (0..=((hi - lo) / step).floor() as usize)
        .map(|i| lo + i as f64 * step)
        .collect();
    radii.push(hi);
    if let Some(path) = o.get("--rc-csv") {
        if Path::new(path) == out {
            bail!("basis output must differ from input curve");
        }
        for row in crate::RotationCurve::from_csv(path)?.rows {
            if (lo..=hi).contains(&row.r_kpc) {
                radii.push(row.r_kpc);
            }
        }
    }
    radii.sort_by(f64::total_cmp);
    radii.dedup();
    let mass = number(&o, "--ref-mass", "1e11")?;
    let mut configs = Vec::new();
    let mut seen = BTreeSet::new();
    for state in value(&o, "--states", DEFAULT_STATES).split(',') {
        for scale in value(&o, "--scales", DEFAULT_SCALES).split(',') {
            let cfg = state_config(state.trim(), scale.trim().parse()?, mass)?;
            if !seen.insert((cfg.n, cfg.l, cfg.m, cfg.a0_star.to_bits())) {
                bail!("duplicate basis state/scale");
            }
            configs.push((state.trim().to_string(), cfg));
        }
    }
    let mut file = output_file(&out)?;
    {
        let mut writer = csv::Writer::from_writer(file.as_file_mut());
        for (id, (state, cfg)) in configs.iter().enumerate() {
            eprintln!(
                "signed basis {}/{}: {} a0={} ({backend})",
                id + 1,
                configs.len(),
                state,
                cfg.a0_star
            );
            let curve = calculate(cfg, &radii, &opts)?;
            if curve.len() != radii.len() || curve.iter().any(|k| !k.is_finite()) {
                bail!("backend returned an invalid signed force curve");
            }
            for (&r, &k) in radii.iter().zip(&curve) {
                writer.serialize(BasisRow {
                    schema_version: 1,
                    basis_id: id,
                    state,
                    n: cfg.n,
                    l: cfg.l,
                    m: cfg.m,
                    a0_star_kpc: cfg.a0_star,
                    ref_mass_msun: mass,
                    gravity_backend: backend,
                    geometry: "azimuthal_mean",
                    normalization: "analytic_full_space",
                    n_radial: opts.nr,
                    n_mu: if backend == "cpu_multipole" {
                        8
                    } else {
                        opts.nmu
                    },
                    n_phi: if backend == "cpu_multipole" {
                        0
                    } else {
                        opts.nphi
                    },
                    softening_kpc: opts.softening,
                    r_kpc: r,
                    v2_signed_kms2: k,
                })?;
            }
        }
        writer.flush()?;
    }
    file.persist(&out).map_err(|e| e.error)?;
    eprintln!("wrote {}", out.display());
    Ok(())
}

#[derive(Deserialize)]
struct WeightedRow {
    n: u32,
    l: u32,
    m: i32,
    a0_star_kpc: f64,
    ref_mass_msun: f64,
    weight: f64,
}

pub fn write_density(args: &[String]) -> Result<()> {
    if args == ["--help"] {
        println!("signed_density --weights FILE --out FILE [--extent 40 --bins 201]");
        return Ok(());
    }
    let o = options(args, &["--weights", "--out", "--extent", "--bins"])?;
    let input = o.get("--weights").context("--weights is required")?;
    let out = Path::new(o.get("--out").context("--out is required")?);
    if Path::new(input) == out {
        bail!("density output must differ from input weights");
    }
    let extent = number(&o, "--extent", "40")?;
    let bins: usize = value(&o, "--bins", "201").parse()?;
    if extent <= 0.0 || !(3..=1001).contains(&bins) {
        bail!("extent must be positive and bins 3..1001");
    }
    let mut reader = csv::Reader::from_path(input)?;
    let mut parts = Vec::new();
    for row in reader.deserialize::<WeightedRow>() {
        let row = row?;
        let cfg = state_config(
            &format!("{}:{}:{}", row.n, row.l, row.m),
            row.a0_star_kpc,
            row.ref_mass_msun,
        )?;
        parts.push((row.weight, cfg));
    }
    let field = if parts.is_empty() {
        None
    } else {
        Some(Superposition::new(parts)?)
    };
    let mut file = output_file(out)?;
    {
        let mut writer = csv::Writer::from_writer(file.as_file_mut());
        writer.write_record(["x_kpc", "z_kpc", "rho_signed_msun_per_kpc3"])?;
        for iz in 0..bins {
            for ix in 0..bins {
                let x = -extent + 2.0 * extent * ix as f64 / (bins - 1) as f64;
                let z = -extent + 2.0 * extent * iz as f64 / (bins - 1) as f64;
                let r = x.hypot(z);
                let theta = if r == 0.0 {
                    0.0
                } else {
                    (z / r).clamp(-1.0, 1.0).acos()
                };
                let rho = match &field {
                    Some(f) => f.rho(r, theta, if x < 0.0 { std::f64::consts::PI } else { 0.0 })?,
                    None => 0.0,
                };
                writer.serialize((x, z, rho))?;
            }
        }
        writer.flush()?;
    }
    file.persist(out).map_err(|e| e.error)?;
    eprintln!("wrote {}", out.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn export_preserves_outward_forces_and_failure_preserves_previous_file() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("basis.csv");
        let args = vec![
            "--out".into(),
            path.to_string_lossy().into(),
            "--states".into(),
            "3d_xy".into(),
            "--scales".into(),
            "1".into(),
            "--r-min".into(),
            "0.5".into(),
            "--r-max".into(),
            "1.5".into(),
            "--r-step".into(),
            "0.5".into(),
        ];
        write_library(&args, "cpu_multipole", |cfg, radii, o| {
            radii
                .iter()
                .map(|r| Ok(orbital_basis::signed_gravity::disk_plane_v2(cfg, *r, o.nr)?))
                .collect()
        })
        .unwrap();
        let before = std::fs::read(&path).unwrap();
        let mut reader = csv::Reader::from_reader(before.as_slice());
        let col = reader
            .headers()
            .unwrap()
            .iter()
            .position(|x| x == "v2_signed_kms2")
            .unwrap();
        let first = reader.records().next().unwrap().unwrap();
        assert!(first[col].parse::<f64>().unwrap() < 0.0);
        assert!(write_library(&args, "cpu_multipole", |_, _, _| bail!("backend failed")).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), before);
    }

    #[test]
    fn density_export_keeps_negative_values_and_accepts_null_fit() {
        let tmp = tempfile::tempdir().unwrap();
        let weights = tmp.path().join("weights.csv");
        let out = tmp.path().join("rho.csv");
        let header = "n,l,m,a0_star_kpc,ref_mass_msun,weight\n";
        std::fs::write(&weights, format!("{header}1,0,0,1,1e11,-1\n")).unwrap();
        let args = vec![
            "--weights".into(),
            weights.to_string_lossy().into(),
            "--out".into(),
            out.to_string_lossy().into(),
            "--extent".into(),
            "1".into(),
            "--bins".into(),
            "3".into(),
        ];
        write_density(&args).unwrap();
        let data = std::fs::read_to_string(&out).unwrap();
        assert!(csv::Reader::from_reader(data.as_bytes())
            .records()
            .all(|r| r.unwrap()[2].parse::<f64>().unwrap() < 0.0));
        std::fs::write(weights, header).unwrap();
        write_density(&args).unwrap();
        let data = std::fs::read_to_string(out).unwrap();
        assert!(csv::Reader::from_reader(data.as_bytes())
            .records()
            .all(|r| r.unwrap()[2].parse::<f64>().unwrap() == 0.0));
    }
}
