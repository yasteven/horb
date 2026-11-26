//! Rebuild the README figures, optionally repeating MW/M31/AGC searches.
use crate::runner::{python, strings, Runner};
use anyhow::{bail, Result};

pub fn run(r: &mut Runner, args: &[&str]) -> Result<()> {
    if args == ["--help"] {
        println!(
            "./run.sh [--release] do readme-highlights [--refit]\n\
Default: rebuild MW/M31/AGC 114905 seven-model comparisons and both density maps,\n\
plus J2141 rotation and lensing comparisons, from recorded fits.\n\
Inputs and scores are checked against the recorded selected weights.\n\
--refit: repeat all 31,931 squared-density and 6,196 signed-psi subsets per galaxy (up to four signed terms).\n\
NFW, Burkert and Einasto are refitted for MW/M31/AGC; QI has no fitted parameters.\n\
J2141 keeps its recorded fits. The single orbital is omitted from these plots.\n\
Outputs: images/*_highlight.png, images/swells/*.png, plots/*.pdf, data/highlights/*.\n\
CPU only; --refit can take several minutes. No Git operations are performed."
        );
        return Ok(());
    }
    if !args.is_empty() && args != ["--refit"] {
        bail!("usage: ./run.sh do readme-highlights [--refit]");
    }
    r.exec(
        python(),
        &strings(&["-c", "import numpy, scipy, matplotlib"]),
        None,
    )?;
    crate::signed::m31_data(r)?;
    // Bound BLAS threading for these small matrices; no global environment changes.
    r.env.insert("OPENBLAS_NUM_THREADS".into(), "1".into());
    r.env.insert("OMP_NUM_THREADS".into(), "1".into());
    for (name, rc, baryons, minimum, maximum, basis_maximum) in [
        (
            "mw",
            "data/milky_way/sofue_2020_standard_rc.csv",
            "data/milky_way/sofue_literature_baryons.csv",
            "5",
            "80",
            "100",
        ),
        (
            "m31",
            "reports/m31_inputs/rotation_curve.csv",
            "reports/m31_inputs/baryons.csv",
            "8.5",
            "36.5",
            "60",
        ),
    ] {
        let basis = format!("reports/{name}_showcase_basis.csv");
        crate::signed::basis(
            r,
            &[
                "signed",
                "cpu",
                "--out",
                &basis,
                "--rc-csv",
                rc,
                "--states",
                "1s,2s,2p_z,3s,3p_z,3d_z2",
                "--scales",
                "0.5,1,2,4,8",
                "--r-min",
                "0.5",
                "--r-max",
                basis_maximum,
                "--r-step",
                "0.5",
                "--n-radial",
                "128",
                "--ref-mass",
                "1e11",
            ],
        )?;
        let prefix = format!("reports/{name}_showcase");
        let (report, weights) = if args == ["--refit"] {
            let mut fit = vec![
                "signed",
                "--basis-csv",
                &basis,
                "--rc-csv",
                rc,
                "--baryons-csv",
                baryons,
                "--out-prefix",
                &prefix,
                "--r-min",
                minimum,
                "--r-max",
                maximum,
                "--max-components",
                "4",
                "--search",
                "exhaustive",
                "--budget",
                "40000",
                "--criterion",
                "chi2",
                "--ridge",
                "0",
                "--max-abs-weight",
                "20",
            ];
            if name == "m31" {
                fit.extend(["--dataset-meta", "reports/m31_inputs/dataset.json"]);
            }
            crate::signed::fit(r, &fit)?;
            (format!("{prefix}.json"), format!("{prefix}_weights.csv"))
        } else {
            (
                format!("data/highlights/{name}_fit.json"),
                format!("data/highlights/{name}_weights.csv"),
            )
        };
        let density = format!("curves/{name}_showcase_density_xz.csv");
        r.curve(
            &[
                "signed_density",
                "--weights",
                &weights,
                "--out",
                &density,
                "--extent",
                "40",
                "--bins",
                "401",
            ],
            None,
        )?;
        let psi_basis = format!("reports/{name}_psi_showcase_basis.csv");
        crate::signed::basis(
            r,
            &[
                "signed",
                "cpu",
                "--density-model",
                "psi",
                "--out",
                &psi_basis,
                "--rc-csv",
                rc,
                "--scales",
                "0.5,1,2,4,8",
                "--r-min",
                "0.5",
                "--r-max",
                basis_maximum,
                "--r-step",
                "0.5",
                "--n-radial",
                "128",
                "--ref-mass",
                "1e11",
            ],
        )?;
        let psi_prefix = format!("reports/{name}_psi_showcase");
        let (psi_report, psi_weights) = if args == ["--refit"]
            || !r
                .root
                .join(format!("data/highlights/{name}_psi_fit.json"))
                .exists()
        {
            crate::signed::fit(
                r,
                &[
                    "signed",
                    "--basis-csv",
                    &psi_basis,
                    "--rc-csv",
                    rc,
                    "--baryons-csv",
                    baryons,
                    "--out-prefix",
                    &psi_prefix,
                    "--r-min",
                    minimum,
                    "--r-max",
                    maximum,
                    "--max-components",
                    "4",
                    "--search",
                    "exhaustive",
                    "--budget",
                    "40000",
                    "--criterion",
                    "chi2",
                    "--ridge",
                    "0",
                    "--max-abs-weight",
                    "20",
                ],
            )?;
            (
                format!("{psi_prefix}.json"),
                format!("{psi_prefix}_weights.csv"),
            )
        } else {
            (
                format!("data/highlights/{name}_psi_fit.json"),
                format!("data/highlights/{name}_psi_weights.csv"),
            )
        };
        let psi_density = format!("curves/{name}_psi_showcase_density_xz.csv");
        r.curve(
            &[
                "signed_density",
                "--weights",
                &psi_weights,
                "--out",
                &psi_density,
                "--extent",
                "40",
                "--bins",
                "401",
            ],
            None,
        )?;
        r.python(
            "make_readme_highlights",
            &[
                "--galaxy",
                name,
                "--fit-report",
                &report,
                "--weights",
                &weights,
                "--basis-csv",
                &basis,
                "--density-csv",
                &density,
                "--psi-fit-report",
                &psi_report,
                "--psi-weights",
                &psi_weights,
                "--psi-basis-csv",
                &psi_basis,
                "--psi-density-csv",
                &psi_density,
            ],
            None,
        )?;
    }
    r.python("make_swells_highlight", &[], None)?;
    crate::agc::run(r, args)?;
    Ok(())
}
