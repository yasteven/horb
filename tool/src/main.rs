mod dump;
mod highlights;
mod runner;
mod signed;
mod workflows;

use anyhow::{bail, Context, Result};
use runner::{cargo, python, strings, Runner};
use std::fs;
use std::io::{self, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const HELP: &str = "HORB · project tool
Usage: ./run.sh [--dry-run] [--release] COMMAND [ARGS...]

  (no arguments)             Interactive menu (help when stdin is not a terminal)
  build [cpu|all|cuda|tool]   Build Rust; all also builds CUDA (default: cpu)
  check [cpu|all|tool|plots] Tests; plots checks Python comparisons (default: cpu)
  fmt                       Format the tool crate
  doctor                    Check Rust, Python plotting packages, CUDA and viewer
  run list                  List project binaries
  run curve [ARGS...]        Run curve_fitter; no args generates the default curve
  run cuda BIN [ARGS...]     Build CUDA library and run one CUDA binary
  run nbody [ARGS...]        Run the existing integrator entry point
  basis signed [cpu|cuda] OPTIONS  Generate signed orbital force curves
  fit signed OPTIONS        Fit orbital choices/counts and signed density weights
  do list                   List reproducible build/run/plot workflows
  do NAME [ARGS...]          Run a workflow (old script stems also accepted)
  data mw                   Regenerate Sofue MW data and standard RC CSV
  data m31                  Prepare Corbelli M31 data and fixed Geehan baryons
  python list               List Python plotting / analysis programs
  python NAME [ARGS...]      Run an existing Python program by stem
  plot NAME [ARGS...]        Run plot_NAME.py (or a full script stem)
  graphs                    List plots/images, numbered newest first
  view [latest|NUMBER|PATH]  Open a graph (default: latest)
  dump [--data] [-o PATH]    Dump current source, including local edits/new files
  help                      Show this help

Aliases: b=build, c=check, r=run, d=dump, source=dump.
Global options go before the command. Put -- after run/plot/python target if desired.
All relative paths resolve from the repository root. No shell command evaluation.
CUDA defaults to sm_87 (Jetson Orin); override with HORB_CUDA_ARCH.
HORB_NVCC, HORB_PYTHON and HORB_VIEWER accept a single executable path.
The dump gets a unique filename under dumps/ and never overwrites an existing file.
";

fn main() {
    if let Err(error) = entry() {
        eprintln!("HORB: {error:#}");
        std::process::exit(1);
    }
}

fn entry() -> Result<()> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .context("missing project root")?
        .canonicalize()?;
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut dry_run = false;
    let mut release = false;
    while let Some(flag) = args.first() {
        match flag.as_str() {
            "--dry-run" => dry_run = true,
            "--release" => release = true,
            _ => break,
        }
        args.remove(0);
    }
    let mut runner = Runner::new(root, dry_run, release);
    if args.is_empty() && io::stdin().is_terminal() {
        return menu(&mut runner);
    }
    dispatch(&mut runner, &args)
}

fn dispatch(r: &mut Runner, args: &[String]) -> Result<()> {
    // Rebuild for each menu action, while nested workflow steps share one build.
    r.reset_cuda();
    let Some(cmd) = args.first().map(String::as_str) else {
        print!("{HELP}");
        return Ok(());
    };
    let rest: Vec<&str> = args[1..].iter().map(String::as_str).collect();
    match cmd {
        "help" | "--help" | "-h" => print!("{HELP}"),
        "build" | "b" => build(r, &rest)?,
        "check" | "c" => check(r, &rest)?,
        "fmt" => {
            no_args(&rest)?;
            r.cargo(&["fmt", "-p", "horb-tool"])?;
        }
        "doctor" => {
            no_args(&rest)?;
            doctor(r)?;
        }
        "do" | "workflow" => {
            let name = rest.first().copied().unwrap_or("list");
            if name == "list" {
                workflows::list();
            } else {
                workflows::run(r, name, &rest[1..])?;
            }
        }
        "data" => match rest.as_slice() {
            ["mw"] => workflows::mw_data(r)?,
            ["m31"] | ["andromeda"] => signed::m31_data(r)?,
            _ => bail!("usage: ./run.sh data mw|m31"),
        },
        "run" | "r" => run(r, &rest)?,
        "basis" => signed::basis(r, &rest)?,
        "fit" => signed::fit(r, &rest)?,
        "python" | "plot" => {
            let name = rest.first().copied().unwrap_or("list");
            if name == "list" {
                list_python(r)?;
            } else {
                let candidate = format!("plot_{}", name.replace('-', "_"));
                let script = if cmd == "plot" && r.python_path(&candidate).is_ok() {
                    &candidate
                } else {
                    name
                };
                if script.trim_end_matches(".py").replace('-', "_") == "plot_cuda_basis_bests_vs_mw"
                {
                    r.cuda_build()?;
                }
                r.dirs()?;
                r.python(script, passthrough(&rest[1..]), None)?;
            }
        }
        "graphs" => {
            no_args(&rest)?;
            graphs(r, true)?;
        }
        "view" => {
            if rest.len() > 1 {
                bail!("usage: ./run.sh view [latest|NUMBER|PATH]");
            }
            view(r, rest.first().copied().unwrap_or("latest"))?;
        }
        "dump" | "d" | "source" => {
            let mut data = false;
            let mut out = None;
            let mut i = 0;
            while i < rest.len() {
                match rest[i] {
                    "--data" => data = true,
                    "-o" | "--output" => {
                        i += 1;
                        out = Some(Path::new(rest.get(i).context("missing output path")?));
                    }
                    "--help" | "-h" => {
                        println!("./run.sh dump [--data] [-o PATH]\nIncludes current Rust/CUDA/Python/shell source, manifests, lockfiles and docs.\n--data also includes input tables under data/. Generated outputs are excluded.");
                        return Ok(());
                    }
                    other => bail!("unknown dump argument '{other}'"),
                }
                i += 1;
            }
            let path = dump::dump(&r.root, out, data, r.dry_run)?;
            println!("{}", path.display());
        }
        _ => bail!("unknown command '{cmd}'; use './run.sh help'"),
    }
    Ok(())
}

fn no_args(args: &[&str]) -> Result<()> {
    if !args.is_empty() {
        bail!("unexpected arguments: {args:?}");
    }
    Ok(())
}
fn passthrough<'a>(args: &'a [&'a str]) -> &'a [&'a str] {
    if args.first() == Some(&"--") {
        &args[1..]
    } else {
        args
    }
}

fn build(r: &mut Runner, args: &[&str]) -> Result<()> {
    if args.len() > 1 {
        bail!("usage: ./run.sh [--release] build [cpu|all|cuda|tool]");
    }
    let target = args.first().copied().unwrap_or("cpu");
    if !["cpu", "all", "cuda", "tool"].contains(&target) {
        bail!("unknown build target '{target}'");
    }
    if target == "all" || target == "cuda" {
        r.cuda_build()?;
    }
    let mut argv = strings(&["build", "--locked"]);
    if r.release {
        argv.push("--release".into());
    }
    match target {
        "all" => argv.push("--workspace".into()),
        "cuda" => argv.extend(strings(&["-p", "cuda_kernels"])),
        "tool" => argv.extend(strings(&["-p", "horb-tool"])),
        _ => argv.extend(strings(&["--workspace", "--exclude", "cuda_kernels"])),
    }
    r.exec(cargo(), &argv, None)
}

fn check(r: &mut Runner, args: &[&str]) -> Result<()> {
    if args.len() > 1 {
        bail!("usage: ./run.sh check [cpu|all|tool|plots]");
    }
    let target = args.first().copied().unwrap_or("cpu");
    if target == "plots" {
        return r.exec(
            python(),
            &strings(&[
                "-m",
                "unittest",
                "discover",
                "-s",
                "scripts/tests",
                "-p",
                "test_*.py",
            ]),
            None,
        );
    }
    if !["cpu", "all", "tool"].contains(&target) {
        bail!("unknown check target '{target}'");
    }
    r.cargo(&["fmt", "-p", "horb-tool", "--", "--check"])?;
    if target == "all" {
        r.cuda_build()?;
    }
    let mut argv = strings(&["test", "--locked"]);
    if r.release {
        argv.push("--release".into());
    }
    match target {
        "all" => argv.push("--workspace".into()),
        "tool" => argv.extend(strings(&["-p", "horb-tool"])),
        _ => argv.extend(strings(&["--workspace", "--exclude", "cuda_kernels"])),
    }
    r.exec(cargo(), &argv, None)
}

fn run(r: &mut Runner, args: &[&str]) -> Result<()> {
    let name = args.first().copied().unwrap_or("list");
    if name == "list" {
        println!("curve (curve_fitter)\nnbody (nbody_integrator; existing placeholder)");
        let mut files = fs::read_dir(r.root.join("cuda_kernels/src/bin"))?
            .collect::<std::io::Result<Vec<_>>>()?;
        files.sort_by_key(|p| p.file_name());
        for f in files {
            if f.path().extension().is_some_and(|e| e == "rs") {
                println!("cuda {}", f.path().file_stem().unwrap().to_string_lossy());
            }
        }
        return Ok(());
    }
    let rest = passthrough(&args[1..]);
    match name {
        "curve" | "curve_fitter" => r.curve(
            if rest.is_empty() {
                &["curve", "3d_z2", "1.5", "1e11"]
            } else {
                rest
            },
            None,
        ),
        "nbody" | "nbody_integrator" => {
            r.exec(cargo(), &r.run_args("nbody_integrator", None, rest), None)
        }
        "cuda" => {
            let bin = rest
                .first()
                .context("usage: ./run.sh run cuda BIN [ARGS...]")?;
            r.cuda(bin, passthrough(&rest[1..]), None)
        }
        _ => bail!("unknown run target '{name}'; use './run.sh run list'"),
    }
}

fn list_python(r: &Runner) -> Result<()> {
    let mut names = Vec::new();
    for dir in ["scripts", "scripts/archive"] {
        if !r.root.join(dir).is_dir() {
            continue;
        }
        for e in fs::read_dir(r.root.join(dir))? {
            let p = e?.path();
            if p.extension().is_some_and(|x| x == "py") {
                names.push(p.file_stem().unwrap().to_string_lossy().into_owned());
            }
        }
    }
    names.sort();
    for name in names {
        println!("{name}");
    }
    Ok(())
}

fn graphs(r: &Runner, show: bool) -> Result<Vec<PathBuf>> {
    let mut found = Vec::new();
    fn collect(dir: &Path, found: &mut Vec<PathBuf>) -> Result<()> {
        if !dir.is_dir() {
            return Ok(());
        }
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            if entry.file_type()?.is_symlink() {
                continue;
            }
            let path = entry.path();
            if path.is_dir() {
                collect(&path, found)?;
            } else if path
                .extension()
                .and_then(|x| x.to_str())
                .is_some_and(|e| ["png", "svg", "pdf", "jpg", "jpeg", "webp"].contains(&e))
            {
                found.push(path);
            }
        }
        Ok(())
    }
    for dir in ["plots", "images"] {
        collect(&r.root.join(dir), &mut found)?;
    }
    found.sort_by(|a, b| {
        let modified = |p: &PathBuf| fs::metadata(p).and_then(|m| m.modified()).ok();
        modified(b).cmp(&modified(a)).then_with(|| a.cmp(b))
    });
    if show {
        if found.is_empty() {
            println!("No graphs yet. Try './run.sh do orbital-views'.");
        }
        for (i, path) in found.iter().enumerate() {
            println!("{:>3}  {}", i + 1, path.strip_prefix(&r.root)?.display());
        }
    }
    Ok(found)
}

fn view(r: &Runner, name: &str) -> Result<()> {
    let path = if name == "latest" || name.parse::<usize>().is_ok() {
        let index = if name == "latest" {
            1
        } else {
            name.parse::<usize>()?
        };
        graphs(r, false)?
            .get(index.checked_sub(1).context("graph numbers start at 1")?)
            .context("graph not found; use './run.sh graphs'")?
            .clone()
    } else {
        r.root.join(name)
    };
    if !path.is_file() {
        bail!("graph does not exist: {}", path.display());
    }
    let viewer = std::env::var_os("HORB_VIEWER").unwrap_or_else(|| {
        if cfg!(target_os = "macos") {
            "open".into()
        } else {
            "xdg-open".into()
        }
    });
    r.exec(viewer, &[path.canonicalize()?.to_string_lossy().into_owned()], None)
        .context("could not open graph; on a headless host use its printed path or set HORB_VIEWER to an executable")
}

fn doctor(r: &Runner) -> Result<()> {
    let checks = [
        ("Rust / Cargo", cargo(), strings(&["--version"])),
        ("Python plotting/fitting", python(), strings(&["-c", "import numpy, scipy, matplotlib, PIL; print('NumPy', numpy.__version__, 'SciPy', scipy.__version__, 'Matplotlib', matplotlib.__version__, 'Pillow', PIL.__version__)"])),
        ("CUDA (optional)", std::env::var_os("HORB_NVCC").unwrap_or_else(|| if Path::new("/usr/local/cuda/bin/nvcc").is_file() { "/usr/local/cuda/bin/nvcc".into() } else { "nvcc".into() }), strings(&["--version"])),
    ];
    for (label, program, args) in checks {
        if r.dry_run {
            eprintln!("check {label}: {}", program.to_string_lossy());
            continue;
        }
        match Command::new(program)
            .args(args)
            .current_dir(&r.root)
            .output()
        {
            Ok(o) if o.status.success() => {
                println!("{label}: OK\n{}", String::from_utf8_lossy(&o.stdout).trim())
            }
            Ok(o) => println!(
                "{label}: unavailable\n{}",
                String::from_utf8_lossy(&o.stderr).trim()
            ),
            Err(e) => println!("{label}: unavailable ({e})"),
        }
    }
    println!(
        "CUDA architecture: {}\nGraph viewer: {}\nProject: {}",
        std::env::var("HORB_CUDA_ARCH").unwrap_or_else(|_| "sm_87".into()),
        std::env::var("HORB_VIEWER").unwrap_or_else(|_| if cfg!(target_os = "macos") {
            "open".into()
        } else {
            "xdg-open".into()
        }),
        r.root.display()
    );
    Ok(())
}

fn menu(r: &mut Runner) -> Result<()> {
    loop {
        println!("\nHORB · Operations\n  1  Build CPU\n  2  Build all (CUDA)\n  3  Check CPU\n  4  Run CUDA Milky Way candidate\n  5  Make orbital density graphs\n  6  List graphs\n  7  View latest graph\n  8  Dump current source\n  9  Dump source + input data\n  s  Fit signed orbital densities (CPU multipoles)\n  g  Fit signed orbital densities (CUDA quadrature)\n  a  Fit Andromeda / M31 (CPU multipoles)\n  d  Dependency checks\n  l  List workflows\n  h  Help\n  q  Quit");
        print!("> ");
        io::stdout().flush()?;
        let mut input = String::new();
        if io::stdin().read_line(&mut input)? == 0 {
            return Ok(());
        }
        let args = match input.trim() {
            "q" | "" => return Ok(()),
            "1" => strings(&["build", "cpu"]),
            "2" => strings(&["build", "all"]),
            "3" => strings(&["check"]),
            "4" => strings(&["do", "cuda-candidate"]),
            "5" => strings(&["do", "orbital-views"]),
            "6" => strings(&["graphs"]),
            "7" => strings(&["view"]),
            "8" => strings(&["dump"]),
            "9" => strings(&["dump", "--data"]),
            "s" => strings(&["do", "signed-fit", "cpu"]),
            "g" => strings(&["do", "signed-fit", "cuda"]),
            "a" => strings(&["do", "m31-signed-fit", "cpu"]),
            "d" => strings(&["doctor"]),
            "l" => strings(&["do", "list"]),
            "h" => strings(&["help"]),
            _ => {
                eprintln!("Choose a menu entry; use './run.sh help' for parameterized commands.");
                continue;
            }
        };
        if let Err(e) = dispatch(r, &args) {
            eprintln!("HORB: {e:#}");
        }
    }
}
