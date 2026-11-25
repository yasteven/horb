use std::process::Command;

#[test]
fn swells_routes_rotation_fit_and_lensing_without_cuda() {
    let output = tool(&[
        "--dry-run",
        "do",
        "swells-fit",
        "--max-components",
        "3",
        "--no-plots",
    ]);
    assert!(output.status.success());
    let plan = String::from_utf8_lossy(&output.stderr);
    assert!(plan.contains("fit_swells.py --max-components 3 --no-plots"));
    assert!(!plan.contains("nvcc"));
    let data = tool(&["--dry-run", "data", "swells"]);
    assert!(data.status.success());
    assert!(String::from_utf8_lossy(&data.stderr).contains("--validate-only"));
    assert!(!tool(&["--dry-run", "do", "swells-fit", "--buget", "100"])
        .status
        .success());
    assert!(!tool(&["--dry-run", "do", "swells-fit", "--budget"])
        .status
        .success());
    let comparison = tool(&[
        "--dry-run",
        "do",
        "swells-compare",
        "--positive-dir",
        "a directory",
        "--no-plots",
    ]);
    assert!(comparison.status.success());
    assert!(String::from_utf8_lossy(&comparison.stderr).contains("compare_swells.py"));
    assert!(
        !tool(&["--dry-run", "do", "swells-compare", "--positive-dir"])
            .status
            .success()
    );
}

#[test]
fn sparc_builds_once_and_forwards_batch_options_without_cuda() {
    let output = tool(&[
        "--dry-run",
        "--release",
        "do",
        "sparc-fit",
        "--galaxies",
        "NGC2403,DDO154",
        "--jobs",
        "2",
        "--name",
        "pilot",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let plan = String::from_utf8(output.stderr).unwrap();
    assert_eq!(plan.matches("--message-format=json").count(), 1);
    assert!(plan.contains("--validate-only"));
    assert!(plan.contains("--curve-bin"));
    assert!(plan.contains("target/release/curve_fitter"));
    assert!(plan.contains("--galaxies NGC2403,DDO154 --jobs 2 --name pilot"));
    assert!(!plan.contains("nvcc"));
    assert!(!tool(&["--dry-run", "do", "sparc-fit", "--jbos", "2"])
        .status
        .success());
    assert!(!tool(&["--dry-run", "do", "sparc-fit", "--jobs"])
        .status
        .success());
    let data = tool(&["--dry-run", "data", "sparc"]);
    assert!(data.status.success());
    assert!(String::from_utf8_lossy(&data.stderr).contains("make_sparc_data.py"));
}

fn tool(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_horb-tool"))
        .env_remove("HORB_CUDA_ARCH")
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn no_input_and_unknown_commands_are_noninteractive() {
    let help = tool(&[]);
    assert!(help.status.success());
    assert!(String::from_utf8_lossy(&help.stdout).contains("dump [--data]"));
    assert!(!tool(&["unknown-command"]).status.success());
    assert!(!tool(&[
        "do",
        "cuda-candidate",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8"
    ])
    .status
    .success());
    assert!(!tool(&["python", "../../arbitrary-file"]).status.success());
    assert!(!tool(&["view", "0"]).status.success());
}

#[test]
fn einasto_parameters_reach_comparisons_but_not_cuda_kernel() {
    for args in [
        vec![
            "--dry-run",
            "do",
            "mw-comparison",
            "1.0",
            "2.5e11",
            "2e11",
            "18",
            "0.25",
        ],
        vec![
            "--dry-run",
            "do",
            "cuda-candidate",
            "1.0",
            "2.5e11",
            "96",
            "80",
            "0.25",
            "18",
            "0.25",
        ],
    ] {
        let result = tool(&args);
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let plan = String::from_utf8(result.stderr).unwrap();
        assert!(plan.contains("5.0 15.0 8.0 80.0 18 0.25"));
        assert!(plan.contains("einasto_rs18_alpha0.25"));
        assert!(plan.contains("score_dm_comparison.py"));
        assert!(plan.contains("_score_inner.csv"));
        assert!(plan.contains("_score_common.csv"));
        if args[2] == "cuda-candidate" {
            let kernel = plan
                .lines()
                .find(|l| l.contains("--bin test_horb_disk_plane_curve"))
                .unwrap();
            assert!(kernel.ends_with("-- 1.0 2.5e11 96 80 0.25"));
            assert!(plan.contains("--cuda-csv"));
        }
    }
    for bad in ["0", "-1", "NaN", "inf", "not-a-number"] {
        let result = tool(&[
            "--dry-run",
            "do",
            "cuda-candidate",
            "1.0",
            "2.5e11",
            "96",
            "80",
            "0.25",
            "15",
            bad,
        ]);
        assert!(!result.status.success());
        assert!(!String::from_utf8_lossy(&result.stderr).contains("nvcc"));
    }
    let plots = tool(&["--dry-run", "check", "plots"]);
    assert!(plots.status.success());
    assert!(String::from_utf8_lossy(&plots.stderr).contains("unittest discover"));
}

#[test]
fn every_workflow_and_legacy_alias_builds_a_dry_run_plan() {
    let list = tool(&["do", "list"]);
    assert!(list.status.success());
    let names = String::from_utf8(list.stdout).unwrap();
    let mut count = 0;
    for line in names.lines().take_while(|line| !line.is_empty()) {
        let name = line.split_whitespace().next().unwrap();
        let mut args = vec!["--dry-run", "do", name];
        if name == "summary" {
            args.push("curves/example.csv");
        }
        let plan = tool(&args);
        assert!(
            plan.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&plan.stderr)
        );
        count += 1;
    }
    assert_eq!(count, 39);
    let alias = tool(&["--dry-run", "do", "scan_cuda_mw_basis.sh"]);
    assert!(alias.status.success());
    assert!(String::from_utf8_lossy(&alias.stderr).contains("--bin scan_horb_cuda_mw_basis"));
}

#[test]
fn cuda_plan_keeps_arguments_and_cpu_build_excludes_cuda() {
    let cpu = tool(&["--dry-run", "build"]);
    let text = String::from_utf8_lossy(&cpu.stderr);
    assert!(text.contains("--exclude cuda_kernels"));
    assert!(!text.contains("nvcc"));
    let cuda = tool(&[
        "--dry-run",
        "--release",
        "run",
        "cuda",
        "write_horb_cuda_curve_basis",
        "--",
        "3d_z2",
        "0.9",
        "2.5e11",
        "64",
        "80",
        "0.25",
    ]);
    assert!(cuda.status.success());
    let text = String::from_utf8_lossy(&cuda.stderr);
    assert!(text.contains("-arch=sm_87"));
    assert!(
        text.contains("--release --bin write_horb_cuda_curve_basis -- 3d_z2 0.9 2.5e11 64 80 0.25")
    );
}

#[test]
fn signed_workflow_routes_raw_force_generation_and_fit_options() {
    for backend in ["cpu", "cuda"] {
        let out = tool(&[
            "--dry-run",
            "do",
            "signed-fit",
            backend,
            "--states",
            "1s,3d_z2",
            "--scales",
            "1,2",
            "--search",
            "anneal",
            "--seed",
            "17",
            "--max-components",
            "2",
        ]);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let plan = String::from_utf8(out.stderr).unwrap();
        assert!(plan.contains("--states 1s,3d_z2"));
        assert!(plan.contains("--search anneal"));
        assert!(plan.contains("--seed 17"));
        assert!(plan.contains("signed_density --weights"));
        assert!(plan.contains("plot_signed_density.py"));
        if backend == "cuda" {
            assert!(plan.contains("--bin write_signed_basis_cuda"));
        } else {
            assert!(plan.contains("-- signed_basis"));
            assert!(!plan.contains("nvcc"));
        }
    }
    assert!(!tool(&["--dry-run", "do", "signed-fit", "--mystery", "1"])
        .status
        .success());
}

#[test]
fn andromeda_routes_its_own_data_and_radius_defaults() {
    let data = tool(&["--dry-run", "data", "m31"]);
    assert!(data.status.success());
    assert!(String::from_utf8_lossy(&data.stderr).contains("make_m31_data.py"));
    for name in ["m31-signed-fit", "andromeda"] {
        let out = tool(&["--dry-run", "do", name]);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let plan = String::from_utf8(out.stderr).unwrap();
        assert!(plan.contains("make_m31_data.py"));
        assert!(plan.contains("--r-min 8.5 --r-max 36.5"));
        assert!(plan.contains("--r-min 0.5 --r-max 60"));
        assert!(plan.contains("--dataset-meta reports/m31_inputs/dataset.json"));
        assert!(plan.contains("plots/m31_signed_cpu_rotation.png"));
        assert!(!plan.contains("sofue"));
    }
    let out = tool(&[
        "--dry-run",
        "do",
        "m31-signed-fit",
        "--name",
        "m31_anneal",
        "--r-min",
        "10",
        "--search",
        "anneal",
    ]);
    assert!(out.status.success());
    let plan = String::from_utf8(out.stderr).unwrap();
    assert!(plan.contains("--r-min 10 --r-max 36.5"));
    assert!(plan.contains("--search anneal"));
    assert!(plan.contains("m31_anneal_rotation.png"));
    let custom = tool(&[
        "--dry-run",
        "do",
        "m31-signed-fit",
        "--rc-csv",
        "custom.csv",
        "--baryons-csv",
        "custom_b.csv",
    ]);
    assert!(custom.status.success());
    let plan = String::from_utf8(custom.stderr).unwrap();
    assert!(!plan.contains("--dataset-meta"));
    assert!(!plan.contains("make_m31_data.py"));
}

#[cfg(unix)]
#[test]
fn child_failure_reaches_cli() {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let child = dir.path().join("failing-cargo");
    std::fs::write(
        &child,
        "#!/bin/sh\nprintf 'incomplete output\\n'\nexit 23\n",
    )
    .unwrap();
    std::fs::set_permissions(&child, std::fs::Permissions::from_mode(0o755)).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_horb-tool"))
        .args(["run", "curve"])
        .env("CARGO", &child)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("23"));
}

#[test]
fn highlights_use_recorded_weights_and_refit_only_when_requested() {
    for refit in [false, true] {
        let mut args = vec!["--dry-run", "--release", "do", "readme-highlights"];
        if refit {
            args.push("--refit");
        }
        let out = tool(&args);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let plan = String::from_utf8(out.stderr).unwrap();
        assert_eq!(plan.matches("make_readme_highlights.py").count(), 2);
        assert_eq!(plan.contains("fit_signed_density.py"), refit);
        assert!(!plan.contains("nvcc"));
        for galaxy in ["mw", "m31"] {
            let prefix = if refit {
                format!("reports/{galaxy}_showcase")
            } else {
                format!("data/highlights/{galaxy}")
            };
            assert!(plan.contains(&format!("signed_density --weights {prefix}_weights.csv")));
        }
        if refit {
            assert!(plan.contains("--criterion chi2 --ridge 0"));
            assert!(plan.contains("--search exhaustive --budget 40000"));
        }
    }
    assert!(
        !tool(&["--dry-run", "do", "readme-highlights", "--mystery"])
            .status
            .success()
    );
}
