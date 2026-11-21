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
    assert!(
        !tool(&["do", "cuda-candidate", "1", "2", "3", "4", "5", "6"])
            .status
            .success()
    );
    assert!(!tool(&["python", "../../arbitrary-file"]).status.success());
    assert!(!tool(&["view", "0"]).status.success());
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
