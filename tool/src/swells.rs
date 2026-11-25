//! Spiral strong lenses with gas kinematics; begin with the tabulated J2141 data.
use crate::runner::Runner;
use anyhow::{bail, Result};

pub fn data(r: &Runner, args: &[&str]) -> Result<()> {
    if !args.is_empty() {
        bail!("usage: ./run.sh data swells");
    }
    r.python("fit_swells", &["--validate-only"], None)
}

pub fn run(r: &mut Runner, args: &[&str]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        return r.python("fit_swells", &["--help"], None);
    }
    let values = [
        "--out-dir",
        "--max-components",
        "--beam-width",
        "--budget",
        "--max-abs-weight",
        "--stellar-mass-offset-dex",
    ];
    let flags = [
        "--nonnegative",
        "--exhaustive",
        "--fit-only",
        "--predict-only",
        "--no-plots",
        "--validate-only",
    ];
    let mut i = 0;
    while i < args.len() {
        if flags.contains(&args[i]) {
            i += 1;
        } else if values.contains(&args[i]) {
            if args.get(i + 1).is_none_or(|value| value.starts_with("--")) {
                bail!("missing value for {}", args[i]);
            }
            i += 2;
        } else {
            bail!(
                "unknown SWELLS option '{}'; use './run.sh do swells-fit --help'",
                args[i]
            );
        }
    }
    r.env.insert("OPENBLAS_NUM_THREADS".into(), "1".into());
    r.env.insert("OMP_NUM_THREADS".into(), "1".into());
    r.python("fit_swells", args, None)
}

pub fn compare(r: &mut Runner, args: &[&str]) -> Result<()> {
    let mut i = 0;
    while i < args.len() {
        match args[i] {
            "--help" | "-h" | "--no-plots" => i += 1,
            "--signed-dir" | "--positive-dir" | "--out-dir" => {
                if args.get(i + 1).is_none_or(|value| value.starts_with("--")) {
                    bail!("missing value for {}", args[i]);
                }
                i += 2;
            }
            other => bail!("unknown SWELLS comparison option '{other}'; use './run.sh do swells-compare --help'"),
        }
    }
    r.env.insert("OPENBLAS_NUM_THREADS".into(), "1".into());
    r.env.insert("OMP_NUM_THREADS".into(), "1".into());
    r.python("compare_swells", args, None)
}
