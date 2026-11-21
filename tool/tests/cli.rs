use std::process::Command;

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
    assert_eq!(count, 33);
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
