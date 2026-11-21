//! Signed density workflow: generate forces once, search weights/subsets, plot.
use crate::runner::{python, strings, Runner};
use anyhow::{bail, Result};
use std::collections::BTreeMap;

const RC: &str = "data/milky_way/sofue_2020_standard_rc.csv";
const BARYONS: &str = "data/milky_way/sofue_literature_baryons.csv";
const M31_RC: &str = "reports/m31_inputs/rotation_curve.csv";
const M31_BARYONS: &str = "reports/m31_inputs/baryons.csv";
const M31_META: &str = "reports/m31_inputs/dataset.json";

fn backend<'a>(args: &'a [&'a str]) -> (&'a str, &'a [&'a str]) {
    match args.first().copied() {
        Some("cpu" | "cuda") => (args[0], &args[1..]),
        _ => ("cpu", args),
    }
}

pub fn basis(r: &mut Runner, args: &[&str]) -> Result<()> {
    if args.first() != Some(&"signed") {
        bail!("usage: ./run.sh basis signed [cpu|cuda] --out FILE [OPTIONS]");
    }
    let (backend, rest) = backend(&args[1..]);
    if backend == "cuda" {
        r.cuda("write_signed_basis_cuda", rest, None)
    } else {
        let mut argv = vec!["signed_basis"];
        argv.extend_from_slice(rest);
        r.curve(&argv, None)
    }
}

pub fn fit(r: &Runner, args: &[&str]) -> Result<()> {
    if args.first() != Some(&"signed") {
        bail!("usage: ./run.sh fit signed --basis-csv FILE --rc-csv FILE --baryons-csv FILE --out-prefix PATH [OPTIONS]");
    }
    r.python("fit_signed_density", &args[1..], None)
}

pub fn workflow(r: &mut Runner, args: &[&str]) -> Result<()> {
    workflow_for(r, args, false)
}

pub fn m31_data(r: &Runner) -> Result<()> {
    r.python("make_m31_data", &[], None)
}

pub fn m31_workflow(r: &mut Runner, args: &[&str]) -> Result<()> {
    workflow_for(r, args, true)
}

fn workflow_for(r: &mut Runner, args: &[&str], m31: bool) -> Result<()> {
    let (backend, args) = backend(args);
    if args == ["--help"] {
        if m31 {
            println!("./run.sh [--release] do m31-signed-fit [cpu|cuda] [OPTIONS]\nM31 defaults: name=m31_signed_cpu/cuda, fit=8.5..36.5 kpc, basis=0.5..60 kpc.\nBundled Corbelli (2010) Table 1 with fixed Geehan (2006) baryons.\nOptions below are shared with the Milky Way workflow; M31 defaults override it.");
        }
        println!(
            "./run.sh [--release] do signed-fit [cpu|cuda] [OPTIONS]\n\
  --name NAME                 Output stem (default mw_signed_cpu/cuda)\n\
  --rc-csv FILE --baryons-csv FILE   Target and fixed baryonic components\n\
  --dataset-meta FILE         Optional provenance with matching input hashes\n\
  --states LIST --scales LIST  Candidate orbitals and a0 values\n\
  --r-min 5 --r-max 80        Fit range in kpc\n\
  --max-components 4 --search beam|anneal|exhaustive\n\
  --beam-width 6 --budget 2000 --seed 42\n\
  --max-abs-weight 20 --ridge 0.0001 --criterion bic|chi2\n\
  --weights signed|nonnegative   Nonnegative is a control comparison\n\
  --n-radial 128 --n-mu 32 --n-phi 64 --softening 0\n\
  --basis-r-min 0.5 --basis-r-max 100 --r-step 0.5 --ref-mass 1e11\n\
  --density-extent 40 --density-bins 201\n\
CPU uses unsoftened axisymmetric multipoles; CUDA uses direct-force quadrature."
        );
        return Ok(());
    }
    let basis_keys = [
        "--states",
        "--scales",
        "--n-radial",
        "--n-mu",
        "--n-phi",
        "--softening",
        "--r-step",
        "--ref-mass",
    ];
    let fit_keys = [
        "--r-min",
        "--r-max",
        "--max-components",
        "--search",
        "--beam-width",
        "--budget",
        "--seed",
        "--max-abs-weight",
        "--ridge",
        "--criterion",
        "--weights",
        "--dataset-meta",
    ];
    let extra_keys = [
        "--name",
        "--rc-csv",
        "--baryons-csv",
        "--basis-r-min",
        "--basis-r-max",
        "--density-extent",
        "--density-bins",
    ];
    if args.len() % 2 != 0 {
        bail!("signed-fit options require --name VALUE pairs; use --help");
    }
    let mut o = BTreeMap::new();
    for pair in args.chunks(2) {
        if !(basis_keys.contains(&pair[0])
            || fit_keys.contains(&pair[0])
            || extra_keys.contains(&pair[0]))
            || o.insert(pair[0], pair[1]).is_some()
        {
            bail!("unknown or duplicate signed-fit option {}", pair[0]);
        }
    }
    if m31 {
        o.entry("--r-min").or_insert("8.5");
        o.entry("--r-max").or_insert("36.5");
        o.entry("--basis-r-max").or_insert("60");
        if !o.contains_key("--rc-csv") && !o.contains_key("--baryons-csv") {
            o.entry("--dataset-meta").or_insert(M31_META);
        }
    }
    let default_name = format!("{}_signed_{backend}", if m31 { "m31" } else { "mw" });
    let name = o.get("--name").copied().unwrap_or(&default_name);
    if name.is_empty()
        || !name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
    {
        bail!("name must contain only letters, digits, underscores or hyphens");
    }
    let rc = o
        .get("--rc-csv")
        .copied()
        .unwrap_or(if m31 { M31_RC } else { RC });
    let baryons =
        o.get("--baryons-csv")
            .copied()
            .unwrap_or(if m31 { M31_BARYONS } else { BARYONS });
    r.exec(
        python(),
        &strings(&["-c", "import numpy, scipy, matplotlib"]),
        None,
    )?;
    if m31 && (!o.contains_key("--rc-csv") || !o.contains_key("--baryons-csv")) {
        m31_data(r)?;
    } else if !m31 && !o.contains_key("--rc-csv") {
        crate::workflows::mw_data(r)?;
    }
    let library = format!("reports/{name}_basis.csv");
    let prefix = format!("reports/{name}");
    let mut library_args = vec![
        "signed",
        backend,
        "--out",
        &library,
        "--rc-csv",
        rc,
        "--r-min",
        o.get("--basis-r-min").copied().unwrap_or("0.5"),
        "--r-max",
        o.get("--basis-r-max").copied().unwrap_or("100"),
    ];
    for key in basis_keys {
        if let Some(v) = o.get(key) {
            library_args.extend([key, v]);
        }
    }
    basis(r, &library_args)?;
    let mut fit_args = vec![
        "signed",
        "--basis-csv",
        &library,
        "--rc-csv",
        rc,
        "--baryons-csv",
        baryons,
        "--out-prefix",
        &prefix,
    ];
    for key in fit_keys {
        if let Some(v) = o.get(key) {
            fit_args.extend([key, v]);
        }
    }
    fit(r, &fit_args)?;
    let density = format!("curves/{name}_density_xz.csv");
    r.curve(
        &[
            "signed_density",
            "--weights",
            &format!("{prefix}_weights.csv"),
            "--out",
            &density,
            "--extent",
            o.get("--density-extent").copied().unwrap_or("40"),
            "--bins",
            o.get("--density-bins").copied().unwrap_or("201"),
        ],
        None,
    )?;
    r.python(
        "plot_signed_density",
        &[&density, "-o", &format!("plots/{name}_density.png")],
        None,
    )?;
    r.python(
        "plot_signed_fit",
        &[
            "--fit-prefix",
            &prefix,
            "-o",
            &format!("plots/{name}_rotation.png"),
        ],
        None,
    )
}
