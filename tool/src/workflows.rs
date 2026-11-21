//! Native orchestration replacing the former shell scripts. Model calculations
//! remain in their existing crates and Python programs.
use crate::runner::Runner;
use anyhow::{bail, Context, Result};
use csv::{Reader, StringRecord, Writer};
use std::fs;

const RC: &str = "data/milky_way/sofue_2020_standard_rc.csv";
const MW: &str = "data/milky_way/sofue_2020_unified_rc.csv";
const BARYONS: &str = "data/milky_way/sofue_literature_baryons.csv";
const SCALES: &[&str] = &["1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "2.0"];
const STATES: &str = "1s,2p_z,2p_x,2p_y,3d_z2,3d_x2_y2,3d_xy,3d_xz,3d_yz";

// (command, former script stem, positional arguments). Old names are aliases.
pub const COMMANDS: &[(&str, &str, &str)] = &[
    (
        "m31-signed-fit",
        "andromeda",
        "[cpu|cuda] [OPTIONS; use --help]",
    ),
    (
        "signed-fit",
        "signed_density_workflow",
        "[cpu|cuda] [OPTIONS; use --help]",
    ),
    ("make-curve", "scripts_make_curve", "[a0=1.5 mass=1e11]"),
    ("density-axes", "scan_dz2_density_axes", ""),
    ("fermi-scales", "scan_dz2_fermi_scales", ""),
    (
        "cuda-candidate",
        "run_cuda_mw_3dz2_candidate",
        "[a0=1.0 mass=2.5e11 n=96 extent=80 softening=0.25 einasto_rs=15.0 alpha=0.18]",
    ),
    (
        "cuda-scan",
        "scan_cuda_mw_3dz2",
        "[rc baryons rmin=5 rmax=25 n=64 extent=80]",
    ),
    ("cuda-convergence", "scan_cuda_mw_3dz2_convergence", ""),
    (
        "cuda-basis",
        "scan_cuda_mw_basis",
        "[rc baryons rmin=5 rmax=25 n=64 extent=80]",
    ),
    (
        "cuda-wavefunction",
        "scan_cuda_mw_wavefunction",
        "[rc baryons rmin=5 rmax=25 n=64 extent=80]",
    ),
    ("cuda-rebuild-plots", "rebuild_cuda_last_step_plots", ""),
    (
        "mw-comparison",
        "run_mw_literature_baryons_model_comparison",
        "[a0=1.0 mass=2.5e11 baseline_mass=2.5e11 einasto_rs=15.0 alpha=0.18]",
    ),
    ("summary", "summarize_total_curve", "CSV"),
    ("orbital-views", "plot_orbital_views", "[a0=1.5 mass=1e11]"),
    ("current-xz", "plot_current_xz", ""),
    ("fermi-views", "plot_fermi_candidate_views", ""),
    ("morphology-scan", "scan_morphology_scales", "[mass=1e11]"),
    (
        "fermi-plots",
        "regenerate_fermi_morphology_plots",
        "[mass=2.5e11]",
    ),
    ("fermi-extent", "scan_fermi_lobe_extent", "[mass=2.5e11]"),
    (
        "fermi-extent40",
        "scan_fermi_lobe_extent_ext40",
        "[mass=2.5e11 extent=40 n=321]",
    ),
    (
        "core-results",
        "regenerate_core_results",
        "(historical spherical/toy-baryon results)",
    ),
    ("legacy-best", "run_best_candidate", ""),
    (
        "legacy-comparison",
        "run_best_candidate_total_model_comparison",
        "",
    ),
    ("legacy-milestone3", "run_milestone3", "[mass=1e11]"),
    ("legacy-milestone4", "run_milestone4", ""),
    ("legacy-milestone5", "run_milestone5", ""),
    ("legacy-sofue-overlay", "run_mw_sofue2020_overlay", ""),
    ("legacy-starter-overlay", "run_mw_starter_overlay", ""),
    (
        "legacy-disk-scan",
        "scan_horb_disk_vs_sofue2020",
        "[a0=1.5]",
    ),
    (
        "legacy-mass-scan",
        "scan_horb_mass_vs_sofue2020",
        "[a0=1.5]",
    ),
    ("legacy-literature-scan", "scan_mw_literature_baryons", ""),
    (
        "legacy-multi-scan",
        "scan_mw_multi_baryons",
        "[rc baryons rmin=5 rmax=25]",
    ),
    (
        "legacy-standard-scan",
        "scan_standard_rc_horb",
        "[rc rmin=5 rmax=25]",
    ),
    ("legacy-score-current", "score_current_sofue_candidates", ""),
    ("legacy-score-ranges", "score_sofue2020_ranges", ""),
];

pub fn list() {
    for (name, _, args) in COMMANDS {
        println!("{name:<26} {args}");
    }
    println!("\nLegacy workflows reproduce historical assumptions. Old .sh stems are accepted.\nUse './run.sh --dry-run do NAME' to inspect a workflow.");
}

fn defaults<'a>(args: &'a [&'a str], values: &[&'a str]) -> Result<Vec<&'a str>> {
    if args.len() > values.len() {
        bail!("too many workflow arguments; use './run.sh do list'");
    }
    Ok(values
        .iter()
        .enumerate()
        .map(|(i, d)| args.get(i).copied().unwrap_or(d))
        .collect())
}

pub fn run(r: &mut Runner, name: &str, args: &[&str]) -> Result<()> {
    let stem = name.trim_end_matches(".sh");
    let name = COMMANDS
        .iter()
        .find(|(n, old, _)| *n == stem || *old == stem)
        .map(|x| x.0)
        .with_context(|| format!("unknown workflow '{name}'; use './run.sh do list'"))?;
    // Validate argument counts before creating outputs or launching commands.
    let maximum = match name {
        "signed-fit" | "m31-signed-fit" => usize::MAX,
        "make-curve" | "orbital-views" => 2,
        "cuda-candidate" => 7,
        "mw-comparison" => 5,
        "cuda-scan" | "cuda-basis" | "cuda-wavefunction" => 6,
        "fermi-extent40" | "legacy-standard-scan" => 3,
        "legacy-multi-scan" => 4,
        "summary" | "morphology-scan" | "fermi-plots" | "fermi-extent" | "legacy-milestone3"
        | "legacy-disk-scan" | "legacy-mass-scan" => 1,
        _ => 0,
    };
    if args.len() > maximum || name == "summary" && args.is_empty() {
        bail!("invalid arguments for {name}; use './run.sh do list'");
    }
    // Environment changes are scoped to one workflow, including nested workflows.
    let saved = r.env.clone();
    let result = run_inner(r, name, args);
    r.env = saved;
    result
}

fn run_inner(r: &mut Runner, name: &str, args: &[&str]) -> Result<()> {
    r.dirs()?;
    match name {
        "signed-fit" => crate::signed::workflow(r, args)?,
        "m31-signed-fit" => crate::signed::m31_workflow(r, args)?,
        "make-curve" => {
            let p = defaults(args, &["1.5", "1e11"])?;
            r.curve(
                &["curve", "3d_z2", p[0], p[1]],
                Some("curves/dz2_curve.csv"),
            )?;
            summary(r, "curves/dz2_curve.csv")?;
        }
        "density-axes" | "fermi-scales" => {
            let scales: &[&str] = if name == "density-axes" {
                &["0.8", "1.0", "1.2", "1.5", "1.7"]
            } else {
                &["0.8", "1.0", "1.2", "1.5", "1.7", "2.5", "8.5"]
            };
            for a in scales {
                let (mode, prefix) = if name == "density-axes" {
                    ("density", "dz2_density_axes")
                } else {
                    ("curve", "dz2")
                };
                let out = format!("curves/{prefix}_a{a}_m1e11.csv");
                r.curve(&[mode, "3d_z2", a, "1e11"], Some(&out))?;
                summary(r, &out)?;
            }
        }
        "cuda-candidate" => cuda_candidate(r, args)?,
        "cuda-scan" | "cuda-basis" | "cuda-wavefunction" => cuda_scan(r, name, args)?,
        "cuda-convergence" => {
            r.env_default("HORB_A0_LIST", "0.9,1.0,1.1");
            r.env_default("HORB_DM_MASS_LIST", "2e11,2.5e11,3e11");
            r.env_default("HORB_SOFTENING_LIST", "0.1,0.25,0.5");
            for n in ["64", "96", "128"] {
                run(r, "cuda-scan", &[RC, BARYONS, "5", "25", n, "80"])?;
            }
        }
        "cuda-rebuild-plots" => cuda_rebuild(r)?,
        "mw-comparison" => mw_comparison(r, args)?,
        "summary" => summary(r, args[0])?,
        "orbital-views" => {
            let p = defaults(args, &["1.5", "1e11"])?;
            orbital_views(r, p[0], p[1])?;
        }
        "fermi-views" => {
            for a in ["1.2", "1.5", "1.7", "2.0"] {
                orbital_views(r, a, "1e11")?;
            }
        }
        "current-xz" => {
            for a in ["1.5", "1.7"] {
                let tag = format!("dz2_xz_a{a}_m1e11");
                r.curve(
                    &["xz", "3d_z2", a, "1e11"],
                    Some(&format!("curves/{tag}.csv")),
                )?;
                r.python(
                    "plot_xz_density",
                    &[
                        &format!("curves/{tag}.csv"),
                        "-o",
                        &format!("plots/{tag}.png"),
                    ],
                    None,
                )?;
            }
        }
        "morphology-scan" => morphology(r, defaults(args, &["1e11"])?[0])?,
        "fermi-plots" => {
            let mass = defaults(args, &["2.5e11"])?[0];
            for a in SCALES {
                let tag = format!("dz2_xz_a{a}_m{mass}");
                r.curve(
                    &["xz", "3d_z2", a, mass],
                    Some(&format!("curves/{tag}.csv")),
                )?;
                density_plots(
                    r,
                    &tag,
                    &[
                        "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.40", "0.50",
                    ],
                    true,
                )?;
            }
        }
        "fermi-extent" => fermi_extent(r, defaults(args, &["2.5e11"])?[0], None)?,
        "fermi-extent40" => {
            let p = defaults(args, &["2.5e11", "40", "321"])?;
            fermi_extent(r, p[0], Some((p[1], p[2])))?;
        }
        "core-results" => core_results(r)?,
        "legacy-best" => legacy_best(r, false)?,
        "legacy-comparison" => legacy_comparison(r)?,
        "legacy-milestone3" => {
            let mass = defaults(args, &["1e11"])?[0];
            morphology(r, mass)?;
            for mode in ["threshold015", "linear", "log"] {
                r.python(
                    "make_morphology_contact_sheet",
                    &["--mode", mode, "--mass", mass],
                    None,
                )?;
            }
            r.python(
                "rank_fermi_candidates",
                &[
                    &format!("reports/morphology_scale_scan_m{mass}.csv"),
                    "--target-min",
                    "8.0",
                    "--target-max",
                    "10.5",
                ],
                None,
            )?;
        }
        "legacy-milestone4" | "legacy-milestone5" => legacy_milestone(r, name)?,
        "legacy-sofue-overlay"
        | "legacy-starter-overlay"
        | "legacy-disk-scan"
        | "legacy-mass-scan" => legacy_scan(r, name, args)?,
        "legacy-literature-scan" | "legacy-multi-scan" | "legacy-standard-scan" => {
            legacy_fit_scan(r, name, args)?
        }
        "legacy-score-current" | "legacy-score-ranges" => score_ranges(r, name)?,
        _ => unreachable!(),
    }
    Ok(())
}

pub fn mw_data(r: &Runner) -> Result<()> {
    r.dirs()?;
    r.python("make_mw_sofue2020_curve", &[], None)?;
    r.python("convert_sofue2020_to_standard_rc", &[], None)
}

fn cuda_candidate(r: &mut Runner, args: &[&str]) -> Result<()> {
    let p = defaults(args, &["1.0", "2.5e11", "96", "80", "0.25", "15.0", "0.18"])?;
    let einasto = einasto_tag(p[5], p[6])?;
    // Match the writer's Rust formatting exactly (also for small exponents).
    let a: f64 = p[0].parse()?;
    let mass: f64 = p[1].parse()?;
    if !a.is_finite() || a <= 0.0 || !mass.is_finite() || mass <= 0.0 {
        bail!("a0 and mass must be finite and positive");
    }
    r.cuda_build()?;
    mw_data(r)?;
    let cuda = format!(
        "curves/cuda_horb_diskplane_a{a:.3}_m{mass:.3e}_n{}.csv",
        p[2]
    );
    let sph = format!(
        "curves/spherical_horb_compare_a{}_m{}_{einasto}.csv",
        p[0], p[1]
    );
    let tag = format!(
        "mw_cuda_horb_diskplane_total_vs_sofue2020_a{}_m{}_n{}_{einasto}",
        p[0], p[1], p[2]
    );
    let plot = format!("plots/{tag}.png");
    r.cuda("test_horb_disk_plane_curve", &p[..5], None)?;
    r.curve(
        &[
            "compare_dm_fixed_baselines",
            "3d_z2",
            p[0],
            p[1],
            p[1],
            "5.0",
            "15.0",
            "8.0",
            "80.0",
            p[5],
            p[6],
        ],
        Some(&sph),
    )?;
    r.python(
        "plot_cuda_total_vs_mw",
        &[
            "--cuda-csv",
            &cuda,
            "--dm-compare-csv",
            &sph,
            "--rc-csv",
            RC,
            "--baryons-csv",
            BARYONS,
            "-o",
            &plot,
        ],
        None,
    )?;
    comparison_scores(r, &sph, Some(&cuda), &tag)
}

fn cuda_scan(r: &mut Runner, name: &str, args: &[&str]) -> Result<()> {
    let p = defaults(args, &[RC, BARYONS, "5", "25", "64", "80"])?;
    r.cuda_build()?;
    mw_data(r)?;
    let (kind, binary) = match name {
        "cuda-basis" => {
            r.env_default("HORB_STATE_LIST", STATES);
            r.env_default("HORB_A0_LIST", "0.6,0.8,0.9,1.0,1.1,1.2,1.4");
            r.env_default("HORB_DM_MASS_LIST", "1.5e11,2e11,2.5e11,3e11");
            r.env_default("HORB_SOFTENING_LIST", "0.25,0.5");
            ("basis", "scan_horb_cuda_mw_basis")
        }
        "cuda-wavefunction" => {
            r.env_default("HORB_WAVE_STATE_LIST", "3d_z2,3d_xy,3d_xz,3d_yz");
            r.env_default("HORB_A0_LIST", "0.8,0.9,1.0,1.1,1.2,1.4");
            r.env_default("HORB_DM_MASS_LIST", "1.5e11,2e11,2.5e11,3e11");
            r.env_default("HORB_SOFTENING_LIST", "0.5");
            r.env_default("HORB_COEFF_LEVELS", "2");
            ("wavefunction", "scan_horb_cuda_mw_wavefunction")
        }
        _ => {
            r.env_default("HORB_A0_LIST", "0.8,0.9,1.0,1.1,1.2,1.3");
            r.env_default("HORB_DM_MASS_LIST", "1.5e11,2e11,2.5e11,3e11,3.5e11");
            r.env_default("HORB_SOFTENING_LIST", "0.1,0.25,0.5");
            ("3dz2", "scan_horb_cuda_mw")
        }
    };
    let out = format!(
        "reports/cuda_mw_{kind}_scan_n{}_r{}_{}.csv",
        p[4], p[2], p[3]
    );
    r.cuda(binary, &[p[0], p[1], &out, p[2], p[3], p[4], p[5]], None)?;
    if !r.dry_run {
        let t = Table::read(r, &out)?;
        t.top("rms_kms", 20)?;
        t.top("chi2_per_point", 20)?;
    }
    Ok(())
}

fn cuda_rebuild(r: &mut Runner) -> Result<()> {
    r.cuda_build()?;
    mw_data(r)?;
    r.env_default("HORB_STATE_LIST", STATES);
    r.env_default("HORB_A0_LIST", "0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4");
    let library = "reports/cuda_orbital_basis_library_n64.csv";
    r.cuda(
        "write_cuda_orbital_basis_library",
        &[library, "64", "80", "0.5", "1e11"],
        None,
    )?;
    r.python(
        "fit_cuda_orbital_basis_wavelets",
        &[
            "--basis-csv",
            library,
            "--rc-csv",
            RC,
            "--baryons-csv",
            BARYONS,
            "--r-min",
            "5",
            "--r-max",
            "25",
            "--max-components",
            "4",
            "--top-basis",
            "24",
            "--metric",
            "chi2",
            "-o",
            "plots/cuda_orbital_basis_wavelet_fit_n64.png",
            "--summary-out",
            "reports/cuda_orbital_basis_wavelet_fit_n64.csv",
        ],
        None,
    )?;
    for (key, source, default) in [
        (
            "HORB_A0_LIST",
            "HORB_WAVE_A0_LIST",
            "0.8,0.9,1.0,1.1,1.2,1.4",
        ),
        (
            "HORB_DM_MASS_LIST",
            "HORB_WAVE_DM_MASS_LIST",
            "1.5e11,2e11,2.5e11,3e11",
        ),
        ("HORB_SOFTENING_LIST", "HORB_WAVE_SOFTENING_LIST", "0.5"),
    ] {
        r.env.insert(
            key.into(),
            std::env::var(source).unwrap_or_else(|_| default.into()),
        );
    }
    run(r, "cuda-wavefunction", &[])?;
    r.python(
        "plot_cuda_wavefunction_scan",
        &[
            "--scan-csv",
            "reports/cuda_mw_wavefunction_scan_n64_r5_25.csv",
            "--rc-csv",
            RC,
            "--baryons-csv",
            BARYONS,
            "--summary-out",
            "reports/cuda_mw_wavefunction_scan_n64_summary.csv",
            "--bar-out",
            "plots/cuda_mw_wavefunction_scan_n64_top_chi2.png",
            "--curve-out",
            "plots/cuda_mw_wavefunction_scan_n64_best_info.png",
            "--metric",
            "chi2_per_point",
            "--top",
            "20",
        ],
        None,
    )?;
    for image in [
        "cuda_orbital_basis_wavelet_fit_n64",
        "cuda_mw_wavefunction_scan_n64_top_chi2",
        "cuda_mw_wavefunction_scan_n64_best_info",
    ] {
        r.copy(
            &format!("plots/{image}.png"),
            &format!("images/{image}.png"),
        )?;
    }
    Ok(())
}

fn mw_comparison(r: &Runner, args: &[&str]) -> Result<()> {
    let p = defaults(args, &["1.0", "2.5e11", "2.5e11", "15.0", "0.18"])?;
    let einasto = einasto_tag(p[3], p[4])?;
    mw_data(r)?;
    let tag = format!(
        "mw_lit_baryons_dm_compare_a{}_m{}_base{}_{einasto}",
        p[0], p[1], p[2]
    );
    let dm = format!("curves/{tag}.csv");
    r.curve(
        &[
            "compare_dm_fixed_baselines",
            "3d_z2",
            p[0],
            p[1],
            p[2],
            "5.0",
            "15.0",
            "8.0",
            "80.0",
            p[3],
            p[4],
        ],
        Some(&dm),
    )?;
    r.python(
        "plot_total_model_comparison_vs_standard_rc",
        &[
            "--dm-compare-csv",
            &dm,
            "--rc-csv",
            RC,
            "--baryons-csv",
            BARYONS,
            "-o",
            &format!("plots/{tag}.png"),
            "--title",
            "Milky Way: spherical HORB vs fixed halo profiles + literature baryons",
        ],
        None,
    )?;
    for (range, max) in [("inner", "25"), ("all", "95.56")] {
        let out = format!("reports/mw_literature_baryons_horb_score_{range}.csv");
        r.curve(
            &[
                "fit_standard_rc_multi_baryons",
                RC,
                BARYONS,
                "3d_z2",
                p[0],
                p[1],
                "5",
                max,
            ],
            Some(&out),
        )?;
    }
    comparison_scores(r, &dm, None, &tag)
}

fn einasto_tag(scale: &str, alpha: &str) -> Result<String> {
    for (name, value) in [("Einasto r_minus2", scale), ("Einasto alpha", alpha)] {
        let number: f64 = value.parse().with_context(|| format!("invalid {name}"))?;
        if !number.is_finite() || number <= 0.0 {
            bail!("{name} must be finite and positive");
        }
    }
    Ok(format!("einasto_rs{scale}_alpha{alpha}"))
}

fn comparison_scores(r: &Runner, dm: &str, cuda: Option<&str>, tag: &str) -> Result<()> {
    for (range, max) in [("inner", "25"), ("common", "80")] {
        let mut args = vec![
            "--dm-compare-csv",
            dm,
            "--rc-csv",
            RC,
            "--baryons-csv",
            BARYONS,
            "--r-min",
            "5",
            "--r-max",
            max,
        ];
        if let Some(path) = cuda {
            args.extend(["--cuda-csv", path]);
        }
        r.python(
            "score_dm_comparison",
            &args,
            Some(&format!("reports/{tag}_score_{range}.csv")),
        )?;
    }
    Ok(())
}

fn density_plots(r: &Runner, tag: &str, thresholds: &[&str], continuous: bool) -> Result<()> {
    let csv = format!("curves/{tag}.csv");
    if continuous {
        for mode in ["log", "linear"] {
            r.python(
                "plot_xz_density",
                &[
                    &csv,
                    "--mode",
                    mode,
                    "--contours",
                    "-o",
                    &format!("plots/{tag}_{mode}.png"),
                    "--title",
                    &format!("{tag}: {mode} density"),
                ],
                None,
            )?;
        }
    }
    for t in thresholds {
        let suffix = format!("{:.2}", t.parse::<f64>()?).replace('.', "");
        r.python(
            "plot_xz_density",
            &[
                &csv,
                "--mode",
                "threshold",
                "--threshold",
                t,
                "--contours",
                "-o",
                &format!("plots/{tag}_threshold{suffix}.png"),
                "--title",
                &format!("{tag}: density threshold {t}"),
            ],
            None,
        )?;
    }
    Ok(())
}

fn orbital_views(r: &Runner, a: &str, mass: &str) -> Result<()> {
    let tag = format!("dz2_xz_a{a}_m{mass}");
    r.curve(
        &["xz", "3d_z2", a, mass],
        Some(&format!("curves/{tag}.csv")),
    )?;
    density_plots(r, &tag, &["0.15"], true)
}

fn morphology(r: &Runner, mass: &str) -> Result<()> {
    let mut report = Table::new("a0_kpc,density_peak_kpc,curve_vmax_radius_kpc,curve_menc_at_vmax_msun,curve_vmax_kms,rho_max_x_kpc,rho_max_z_kpc,rho_max_r_kpc,rho_max_msun_per_kpc3");
    for a in ["1.0"].iter().chain(SCALES) {
        let curve = format!("curves/dz2_curve_a{a}_m{mass}.csv");
        r.curve(&["curve", "3d_z2", a, mass], Some(&curve))?;
        orbital_views(r, a, mass)?;
        if !r.dry_run {
            let c = Table::read(r, &curve)?;
            let peak = c.maximum("v_circ_kms")?;
            let d = Table::read(r, &format!("curves/dz2_xz_a{a}_m{mass}.csv"))?;
            let rho = d.maximum("rho_Msun_per_kpc3")?;
            let density_peak = format!("{:.6}", 6.0 * a.parse::<f64>()?);
            report.rows.push(StringRecord::from(vec![
                *a,
                &density_peak,
                &peak[0],
                &peak[1],
                &peak[2],
                &rho[d.col("x_kpc")?],
                &rho[d.col("z_kpc")?],
                &rho[d.col("r_kpc")?],
                &rho[d.col("rho_Msun_per_kpc3")?],
            ]));
        }
    }
    report.save(r, &format!("reports/morphology_scale_scan_m{mass}.csv"))
}

fn fermi_extent(r: &mut Runner, mass: &str, grid: Option<(&str, &str)>) -> Result<()> {
    let suffix = grid
        .map(|(ext, _)| format!("_ext{ext}"))
        .unwrap_or_default();
    if let Some((ext, n)) = grid {
        r.env.insert("HORB_XZ_EXTENT_KPC".into(), ext.into());
        r.env.insert("HORB_XZ_N".into(), n.into());
    }
    let mut report = Table::default();
    for a in SCALES {
        let tag = format!("dz2_xz_a{a}_m{mass}{suffix}");
        let csv = format!("curves/{tag}.csv");
        r.curve(&["xz", "3d_z2", a, mass], Some(&csv))?;
        let mut args = vec![
            csv.as_str(),
            "--thresholds",
            "0.05",
            "0.10",
            "0.15",
            "0.20",
            "0.25",
        ];
        if grid.is_some() {
            args.extend(["0.30", "0.35", "0.40", "0.45"]);
        }
        args.push("0.50");
        report.extend(&r.python_capture("measure_xz_lobe_extent", &args)?)?;
        density_plots(
            r,
            &tag,
            if grid.is_some() {
                &["0.15", "0.20", "0.25", "0.30", "0.40", "0.50"]
            } else {
                &["0.15"]
            },
            false,
        )?;
    }
    report.save(
        r,
        &format!("reports/fermi_lobe_extent_scan_m{mass}{suffix}.csv"),
    )
}

fn legacy_best(r: &mut Runner, core: bool) -> Result<()> {
    mw_data(r)?;
    let total = if core {
        "curves/best_candidate_total_a1.5_m2.5e11_disk1e11.csv"
    } else {
        "curves/best_candidate_total.csv"
    };
    r.curve(
        &[
            "total", "3d_z2", "1.5", "2.5e11", "1e11", "3.0", "1e10", "0.7",
        ],
        Some(total),
    )?;
    r.python(
        "plot_total_vs_mw",
        &[
            total,
            MW,
            "-o",
            "plots/best_candidate_total_vs_sofue2020.png",
            "--title",
            "Historical spherical HORB candidate with toy baryons",
        ],
        None,
    )?;
    for (label, max) in [("inner", "25"), ("all", "95.56")] {
        let label = if core {
            format!("{label}_5_{}", if max == "25" { "25" } else { "95" })
        } else {
            label.into()
        };
        r.python(
            "score_total_vs_mw",
            &[total, MW, "--r-min", "5", "--r-max", max],
            Some(&format!("reports/best_candidate_score_{label}.csv")),
        )?;
    }
    let tag = if core {
        "best_candidate_xz_a1.5_m2.5e11_ext40"
    } else {
        "best_candidate_xz_ext40"
    };
    let csv = format!("curves/{tag}.csv");
    r.env.insert("HORB_XZ_EXTENT_KPC".into(), "40".into());
    r.env.insert("HORB_XZ_N".into(), "321".into());
    r.curve(&["xz", "3d_z2", "1.5", "2.5e11"], Some(&csv))?;
    density_plots(r, tag, &["0.15", "0.20", "0.25", "0.50"], true)?;
    r.python(
        "measure_xz_lobe_extent",
        &[
            &csv,
            "--thresholds",
            "0.10",
            "0.15",
            "0.20",
            "0.25",
            "0.30",
            "0.40",
            "0.50",
        ],
        Some("reports/best_candidate_lobe_extent.csv"),
    )
}

fn core_results(r: &mut Runner) -> Result<()> {
    legacy_best(r, true)?;
    r.python(
        "plot_mw_curve",
        &[MW, "-o", "plots/sofue_2020_unified_rc.png"],
        None,
    )?;
    run(r, "fermi-extent40", &[])?;
    for t in ["015", "020", "025", "050"] {
        r.python(
            "make_fermi_ext40_contact_sheet",
            &["--mass", "2.5e11", "--extent", "40", "--threshold", t],
            None,
        )?;
    }
    run(r, "legacy-disk-scan", &["1.5"])?;
    run(r, "legacy-score-current", &[])
}

fn legacy_comparison(r: &Runner) -> Result<()> {
    mw_data(r)?;
    let total = "curves/best_candidate_total_a1.5_m2.5e11_disk1e11.csv";
    let dm = "curves/best_candidate_dm_compare_baselines.csv";
    r.curve(
        &[
            "total", "3d_z2", "1.5", "2.5e11", "1e11", "3.0", "1e10", "0.7",
        ],
        Some(total),
    )?;
    r.curve(
        &[
            "compare_dm_fixed_baselines",
            "3d_z2",
            "1.5",
            "2.5e11",
            "2.5e11",
            "5.0",
            "15.0",
            "8.0",
            "80.0",
        ],
        Some(dm),
    )?;
    let common = ["--total-csv", total, "--dm-compare-csv", dm, "--mw-csv", RC];
    let mut plot = common.to_vec();
    plot.extend(["-o", "plots/best_candidate_total_models_vs_sofue2020.png"]);
    r.python("plot_total_model_comparison_vs_mw", &plot, None)?;
    for (label, max) in [("inner", "25"), ("all", "95.56")] {
        let mut score = common.to_vec();
        score.extend(["--r-min", "5", "--r-max", max]);
        r.python(
            "score_total_model_comparison_vs_mw",
            &score,
            Some(&format!(
                "reports/best_candidate_total_models_score_{label}.csv"
            )),
        )?;
    }
    Ok(())
}

fn legacy_milestone(r: &Runner, name: &str) -> Result<()> {
    for a in ["1.4", "1.5", "1.6", "1.7"] {
        if name == "legacy-milestone4" {
            let tag = format!("dz2_total_a{a}_m1e11_disk6e10_bulge1e10");
            let csv = format!("curves/{tag}.csv");
            r.curve(
                &["total", "3d_z2", a, "1e11", "6e10", "3.0", "1e10", "0.7"],
                Some(&csv),
            )?;
            summary(r, &csv)?;
            r.python(
                "plot_total_curve",
                &[&csv, "-o", &format!("plots/{tag}.png")],
                None,
            )?;
        } else {
            let tag = format!("compare_dm_horb_dz2_a{a}_m1e11_vs_baselines");
            let csv = format!("curves/{tag}.csv");
            r.curve(
                &[
                    "compare_dm",
                    "3d_z2",
                    a,
                    "1e11",
                    "5.0",
                    "15.0",
                    "8.0",
                    "80.0",
                ],
                Some(&csv),
            )?;
            r.python(
                "plot_dm_comparison",
                &[&csv, "-o", &format!("plots/{tag}.png")],
                None,
            )?;
            r.python(
                "summarize_dm_comparison",
                &[&csv],
                Some(&format!("reports/{tag}_summary.csv")),
            )?;
        }
    }
    Ok(())
}

fn legacy_scan(r: &Runner, name: &str, args: &[&str]) -> Result<()> {
    mw_data(r)?;
    let a = args.first().copied().unwrap_or("1.5");
    let (scales, masses, disks, mw, out): (Vec<&str>, Vec<&str>, Vec<&str>, &str, String) =
        match name {
            "legacy-disk-scan" => (
                vec![a],
                vec!["1.5e11", "2e11", "2.5e11", "3e11"],
                vec!["6e10", "8e10", "1e11"],
                MW,
                format!("reports/score_horb_disk_scan_a{a}_vs_sofue2020.csv"),
            ),
            "legacy-mass-scan" => (
                vec![a],
                vec![
                    "5e10", "6e10", "7e10", "8e10", "9e10", "1e11", "1.2e11", "1.5e11",
                ],
                vec!["6e10"],
                MW,
                format!("reports/score_horb_mass_scan_a{a}_vs_sofue2020.csv"),
            ),
            "legacy-starter-overlay" => {
                r.python("make_mw_starter_curve", &[], None)?;
                (
                    vec!["1.4", "1.5", "1.6", "1.7"],
                    vec!["5e10"],
                    vec!["6e10"],
                    "data/milky_way/mw_starter_total_curve.csv",
                    "reports/score_horb_total_vs_mw_starter.csv".into(),
                )
            }
            _ => (
                vec!["1.4", "1.5", "1.6", "1.7"],
                vec!["4e10", "5e10", "6e10", "7e10"],
                vec!["6e10"],
                MW,
                "reports/score_horb_total_vs_sofue2020.csv".into(),
            ),
        };
    let mut table = Table::default();
    for a in scales {
        for mass in &masses {
            for disk in &disks {
                let tag = format!("dz2_total_a{a}_m{mass}_disk{disk}_bulge1e10");
                let csv = format!("curves/{tag}.csv");
                r.curve(
                    &["total", "3d_z2", a, mass, disk, "3.0", "1e10", "0.7"],
                    Some(&csv),
                )?;
                r.python(
                    "plot_total_vs_mw",
                    &[
                        &csv,
                        mw,
                        "-o",
                        &format!(
                            "plots/{tag}_vs_{}.png",
                            if name == "legacy-starter-overlay" {
                                "mw_starter"
                            } else {
                                "sofue2020"
                            }
                        ),
                    ],
                    None,
                )?;
                table.extend(&r.python_capture(
                    "score_total_vs_mw",
                    &[&csv, mw, "--r-min", "5", "--r-max", "25"],
                )?)?;
            }
        }
    }
    table.save(r, &out)?;
    if !r.dry_run {
        table.top("chi2_per_point", 20)?;
    }
    Ok(())
}

fn legacy_fit_scan(r: &Runner, name: &str, args: &[&str]) -> Result<()> {
    mw_data(r)?;
    let mut report = Table::default();
    let out;
    if name == "legacy-multi-scan" {
        let p = defaults(args, &[RC, BARYONS, "5", "25"])?;
        out = "reports/mw_multi_baryons_horb_scan.csv";
        for a in ["1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6"] {
            for m in ["2e11", "2.5e11", "3e11", "3.5e11", "4e11", "4.5e11", "5e11"] {
                report.extend(&r.curve_capture(&[
                    "fit_standard_rc_multi_baryons",
                    p[0],
                    p[1],
                    "3d_z2",
                    a,
                    m,
                    p[2],
                    p[3],
                ])?)?;
            }
        }
    } else {
        let p = defaults(args, &[RC, "5", "25"])?;
        let (scales, masses, disks, ds, bm, bs): (&[&str], &[&str], &[&str], &str, &str, &str) =
            if name == "legacy-literature-scan" {
                out = "reports/mw_literature_baryons_horb_scan.csv";
                (
                    &["1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8"],
                    &[
                        "2e11", "2.5e11", "3e11", "3.5e11", "4e11", "5e11", "6e11", "8e11",
                    ],
                    &["3.41e10"],
                    "3.19",
                    "1.652e10",
                    "0.522",
                )
            } else {
                out = "reports/standard_rc_scan_horb.csv";
                (
                    &["1.3", "1.4", "1.5", "1.6", "1.7"],
                    &["2e11", "2.5e11", "3e11"],
                    &["8e10", "1e11", "1.2e11"],
                    "3.0",
                    "1e10",
                    "0.7",
                )
            };
        for a in scales {
            for m in masses {
                for disk in disks {
                    report.extend(&r.curve_capture(&[
                        "fit_standard_rc",
                        p[0],
                        "3d_z2",
                        a,
                        m,
                        disk,
                        ds,
                        bm,
                        bs,
                        p[1],
                        p[2],
                    ])?)?;
                }
            }
        }
    }
    report.save(r, out)?;
    if !r.dry_run {
        report.top("rms_kms", 20)?;
    }
    Ok(())
}

fn score_ranges(r: &Runner, name: &str) -> Result<()> {
    let current = name == "legacy-score-current";
    let out = if current {
        "reports/score_current_sofue_candidates_ranges.csv"
    } else {
        "reports/score_horb_total_vs_sofue2020_ranges.csv"
    };
    if r.dry_run {
        eprintln!("score matching curves/dz2_total_*.csv by radius range -> {out}");
        return Ok(());
    }
    let mut paths = fs::read_dir(r.root.join("curves"))?.collect::<std::io::Result<Vec<_>>>()?;
    paths.sort_by_key(|p| p.file_name());
    let mut report = Table::new("model_csv,range_name,r_min,r_max,n,rms_kms,chi2,chi2_per_point");
    for entry in paths {
        let file = entry.file_name().to_string_lossy().into_owned();
        if !file.ends_with("_bulge1e10.csv")
            || !file.starts_with(if current {
                "dz2_total_a1.5_m"
            } else {
                "dz2_total_a"
            })
            || !file.contains("_disk")
        {
            continue;
        }
        if !current && !file.contains("_disk6e10_") {
            continue;
        }
        let mut ranges = vec![("inner", "5", "25")];
        if current {
            ranges.push(("mid", "10", "40"));
        }
        ranges.extend([("outer", "25", "95.56"), ("all", "5", "95.56")]);
        for (label, min, max) in ranges {
            let mut score = Table::default();
            score.extend(&r.python_capture(
                "score_total_vs_mw",
                &[
                    &format!("curves/{file}"),
                    MW,
                    "--r-min",
                    min,
                    "--r-max",
                    max,
                ],
            )?)?;
            for row in &score.rows {
                let mut fields = vec![row[0].to_string(), label.into()];
                fields.extend(row.iter().skip(1).map(String::from));
                report.rows.push(StringRecord::from(fields));
            }
        }
    }
    if report.rows.is_empty() {
        bail!("no matching total curves; run a legacy overlay or disk scan first");
    }
    report.save(r, out)
}

fn summary(r: &Runner, path: &str) -> Result<()> {
    if r.dry_run {
        eprintln!("summarize maxima in {path}");
        return Ok(());
    }
    let table = Table::read(r, path)?;
    println!("file: {path}\nrows: {}", table.rows.len());
    for (i, name) in table.header.iter().enumerate().skip(1) {
        if let Ok(row) = table.maximum(name) {
            println!(
                "{name}: max={} at {}={}",
                &row[i], &table.header[0], &row[0]
            );
        }
    }
    Ok(())
}

#[derive(Default)]
struct Table {
    header: StringRecord,
    rows: Vec<StringRecord>,
}
impl Table {
    fn new(header: &str) -> Self {
        Self {
            header: StringRecord::from(header.split(',').collect::<Vec<_>>()),
            rows: Vec::new(),
        }
    }
    fn read(r: &Runner, path: &str) -> Result<Self> {
        let mut t = Self::default();
        t.extend(&fs::read_to_string(r.root.join(path)).with_context(|| format!("read {path}"))?)?;
        Ok(t)
    }
    fn extend(&mut self, text: &str) -> Result<()> {
        if text.trim().is_empty() {
            return Ok(());
        } // Dry runs produce no child output.
        let mut reader = Reader::from_reader(text.as_bytes());
        let header = reader.headers()?.clone();
        if self.header.is_empty() {
            self.header = header;
        } else if self.header != header {
            bail!("CSV headers changed while combining reports");
        }
        for row in reader.records() {
            self.rows.push(row?);
        }
        Ok(())
    }
    fn col(&self, name: &str) -> Result<usize> {
        self.header
            .iter()
            .position(|s| s == name)
            .with_context(|| format!("missing CSV column '{name}'"))
    }
    fn maximum(&self, name: &str) -> Result<&StringRecord> {
        let col = self.col(name)?;
        self.rows
            .iter()
            .filter_map(|row| {
                row[col]
                    .parse::<f64>()
                    .ok()
                    .filter(|v| v.is_finite())
                    .map(|value| (row, value))
            })
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|x| x.0)
            .context("no finite rows")
    }
    fn save(&self, r: &Runner, path: &str) -> Result<()> {
        eprintln!("write report: {path}");
        if r.dry_run {
            return Ok(());
        }
        let mut w = Writer::from_writer(Vec::new());
        w.write_record(&self.header)?;
        for row in &self.rows {
            w.write_record(row)?;
        }
        r.write(path, &String::from_utf8(w.into_inner()?)?)
    }
    fn top(&self, key: &str, count: usize) -> Result<()> {
        let col = self.col(key)?;
        let mut rows: Vec<_> = self
            .rows
            .iter()
            .filter_map(|row| {
                row[col]
                    .parse::<f64>()
                    .ok()
                    .filter(|n| n.is_finite())
                    .map(|score| (score, row))
            })
            .collect();
        rows.sort_by(|a, b| a.0.total_cmp(&b.0));
        println!("Best by {key}:");
        let mut writer = Writer::from_writer(std::io::stdout());
        writer.write_record(&self.header)?;
        for (_, row) in rows.into_iter().take(count) {
            writer.write_record(row)?;
        }
        writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn aggregation_preserves_quoted_fields_and_rejects_changed_headers() {
        let mut t = Table::default();
        t.extend("name,rms\n\"a,b\",2\n").unwrap();
        t.extend("name,rms\nc,1\n").unwrap();
        assert_eq!(&t.rows[0][0], "a,b");
        assert_eq!(&t.maximum("rms").unwrap()[0], "a,b");
        assert!(t.extend("different,rms\nx,1\n").is_err());
    }
}
