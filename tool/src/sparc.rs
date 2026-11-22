//! Build once, then run the offline, resumable per-galaxy SPARC workflow.
use crate::runner::{cargo, strings, Runner};
use anyhow::{bail, Context, Result};
use std::path::PathBuf;

pub fn data(r: &Runner, args: &[&str]) -> Result<()> {
    r.python("make_sparc_data", args, None)
}

pub fn run(r: &mut Runner, args: &[&str]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        return r.python("fit_sparc", &["--help"], None);
    }
    // Reject misspelled switches even when --dry-run does not execute Python.
    let values = [
        "--name",
        "--galaxies",
        "--jobs",
        "--search",
        "--budget",
        "--beam-width",
        "--seed",
        "--max-components",
        "--max-abs-weight",
        "--criterion",
        "--ridge",
        "--weights",
        "--ml-disk",
        "--ml-bulge",
        "--states",
        "--scale-factors",
        "--scales",
        "--ref-mass",
        "--n-radial",
        "--grid-points",
        "--plots",
        "--density-bins",
    ];
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--force" {
            i += 1;
        } else if values.contains(&args[i]) {
            if args.get(i + 1).is_none_or(|value| value.starts_with("--")) {
                bail!("missing value for {}", args[i]);
            }
            i += 2;
        } else {
            bail!(
                "unknown SPARC option '{}'; use './run.sh do sparc-fit --help'",
                args[i]
            );
        }
    }
    let mut validate = args.to_vec();
    validate.push("--validate-only");
    r.python("fit_sparc", &validate, None)?;
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
    let executable = if r.dry_run {
        // A display-only placeholder; real runs use Cargo's executable path,
        // including custom CARGO_TARGET_DIR/target configuration.
        r.root
            .join("target")
            .join(if r.release { "release" } else { "debug" })
            .join("curve_fitter")
    } else {
        cargo_executable(&output)?
    };
    let executable = executable.to_str().context("non-UTF8 curve_fitter path")?;
    let mut batch = args.to_vec();
    batch.extend(["--curve-bin", executable]);
    r.python("fit_sparc", &batch, None)
}

fn cargo_executable(output: &str) -> Result<PathBuf> {
    output
        .lines()
        .filter_map(|line| serde_json::from_str::<serde_json::Value>(line).ok())
        .filter(|v| v["reason"] == "compiler-artifact" && v["target"]["name"] == "curve_fitter")
        .filter_map(|v| v["executable"].as_str().map(PathBuf::from))
        .last()
        .context("Cargo did not report the curve_fitter executable")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_actual_cargo_executable_with_custom_target_directory() {
        let messages = r#"{"reason":"build-finished","success":true}
{"reason":"compiler-artifact","target":{"name":"curve_fitter"},"executable":"/tmp/custom target/release/curve_fitter"}"#;
        assert_eq!(
            cargo_executable(messages).unwrap(),
            PathBuf::from("/tmp/custom target/release/curve_fitter")
        );
        assert!(cargo_executable("{}").is_err());
    }
}
