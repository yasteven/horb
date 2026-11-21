use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_curve_fitter"))
        .args(args)
        .output()
        .unwrap()
}

fn csv(args: &[&str]) -> (Vec<String>, Vec<Vec<f64>>) {
    let out = run(args);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = String::from_utf8(out.stdout).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap().split(',').map(String::from).collect();
    let rows = lines
        .map(|line| line.split(',').map(|x| x.parse().unwrap()).collect())
        .collect();
    (header, rows)
}

#[test]
fn standalone_einasto_has_default_and_configurable_shape() {
    let (header, default) = csv(&["halo", "einasto", "15", "2.5e11", "80"]);
    let (_, explicit) = csv(&["halo", "einasto", "15", "2.5e11", "80", "0.18"]);
    let (_, different) = csv(&["halo", "einasto", "15", "2.5e11", "80", "0.5"]);
    assert_eq!(header.len(), 4); // Preserve the standalone halo schema.
    assert_eq!(default, explicit);
    assert_ne!(default[50], different[50]);
    assert!(default.iter().flatten().all(|x| x.is_finite() && *x > 0.0));
}

#[test]
fn comparison_keeps_existing_columns_and_normalizes_every_baseline() {
    let args = [
        "compare_dm_fixed_baselines",
        "1s",
        "1",
        "1e11",
        "2.5e11",
        "5",
        "15",
        "8",
        "80",
    ];
    let (header, default) = csv(&args);
    assert_eq!(
        &header[..9],
        [
            "r_kpc",
            "v_horb_kms",
            "v_piso_kms",
            "v_nfw_kms",
            "v_burkert_kms",
            "m_horb_msun",
            "m_piso_msun",
            "m_nfw_msun",
            "m_burkert_msun"
        ]
    );
    assert_eq!(
        &header[9..],
        [
            "v_einasto_kms",
            "m_einasto_msun",
            "einasto_r_minus2_kpc",
            "einasto_alpha",
            "einasto_rho_minus2_msun_per_kpc3"
        ]
    );
    let mut custom_args = args.to_vec();
    custom_args.extend(["20", "0.3"]);
    let (_, custom) = csv(&custom_args);
    for (old, new) in default.iter().zip(&custom) {
        assert_eq!(&old[..9], &new[..9]);
        assert_eq!((old[11], old[12]), (15.0, 0.18));
        assert_eq!((new[11], new[12]), (20.0, 0.3));
        assert!(new.iter().all(|x| x.is_finite() && *x > 0.0));
    }
    assert_ne!(default[100][9], custom[100][9]);
    for rows in [&default, &custom] {
        let end = rows.last().unwrap();
        assert_eq!(end[0], 80.0);
        for col in [6, 7, 8, 10] {
            assert!((end[col] / 2.5e11 - 1.0).abs() < 1e-6);
        }
        assert!((end[5] / 1e11 - 1.0).abs() < 1e-6);
    }
}

#[test]
fn original_compare_form_uses_its_mass_and_nfw_scale_for_einasto() {
    let (_, rows) = csv(&["compare_dm", "1s", "1", "1e11", "5", "18", "8", "80"]);
    let end = rows.last().unwrap();
    assert!((end[10] / 1e11 - 1.0).abs() < 1e-6);
    assert_eq!(end[11], 18.0);
    assert_eq!(end[12], 0.18);
}

#[test]
fn invalid_shape_is_rejected_without_partial_csv() {
    for bad in ["0", "-1", "NaN", "inf", "garbage"] {
        for args in [
            vec!["halo", "einasto", "15", "2.5e11", "80", bad],
            vec![
                "compare_dm",
                "1s",
                "1",
                "1e11",
                "5",
                "15",
                "8",
                "80",
                "15",
                bad,
            ],
        ] {
            let result = run(&args);
            assert!(!result.status.success());
            assert!(result.stdout.is_empty());
        }
    }
    assert!(!run(&["halo", "nfw", "15", "2.5e11", "80", "0.18"])
        .status
        .success());
}
