//! AGC 114905: attributed digitized observations, paired HORB fits and baselines.
use crate::runner::{cargo, strings, Runner};
use anyhow::{bail, Context, Result};

pub fn run(r: &Runner, args: &[&str]) -> Result<()> {
    if args == ["--help"] {
        println!("./run.sh --release do agc114905-fit [--refit]\nFit squared-density and signed-psi HORB, QI and classical baselines to the\n2024 Figure 7 digitized data. Include zero orbital coefficients (baryons alone).\nDefault reuses checked recorded HORB fits; --refit repeats both subset searches.");
        return Ok(());
    }
    if !args.is_empty() && args != ["--refit"] {
        bail!("usage: do agc114905-fit [--refit]");
    }
    let mut build = strings(&[
        "build",
        "--locked",
        "-p",
        "curve_fitter",
        "--bin",
        "curve_fitter",
        "--message-format=json",
    ]);
    if r.release {
        build.push("--release".into());
    }
    let output = r.capture(cargo(), &build)?;
    let bin = if r.dry_run {
        r.root
            .join("target")
            .join(if r.release { "release" } else { "debug" })
            .join("curve_fitter")
    } else {
        crate::sparc::cargo_executable(&output)?
    };
    let mut command = vec![
        "--curve-bin",
        bin.to_str().context("non-UTF8 executable path")?,
    ];
    command.extend_from_slice(args);
    r.python("fit_agc114905", &command, None)
}
