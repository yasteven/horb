use anyhow::{bail, Context, Result};
use std::collections::BTreeMap;
use std::ffi::{OsStr, OsString};
use std::fs::{self, OpenOptions};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub struct Runner {
    pub root: PathBuf,
    pub dry_run: bool,
    pub release: bool,
    pub env: BTreeMap<String, String>,
    cuda_ready: bool,
}

pub fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|s| s.to_string()).collect()
}

impl Runner {
    pub fn reset_cuda(&mut self) {
        self.cuda_ready = false;
    }

    pub fn new(root: PathBuf, dry_run: bool, release: bool) -> Self {
        Self {
            root,
            dry_run,
            release,
            env: BTreeMap::new(),
            cuda_ready: false,
        }
    }

    fn command(&self, program: impl AsRef<OsStr>, args: &[String]) -> Result<Command> {
        let mut cmd = Command::new(program);
        cmd.args(args).current_dir(&self.root).envs(&self.env);
        let mut paths = vec![self.root.join("cuda_kernels/cuda")];
        if let Some(value) = std::env::var_os("LD_LIBRARY_PATH") {
            paths.extend(std::env::split_paths(&value).filter(|p| !p.as_os_str().is_empty()));
        }
        cmd.env("LD_LIBRARY_PATH", std::env::join_paths(paths)?);
        // Plot commands save figures and work on a headless SSH host.
        cmd.env("MPLBACKEND", "Agg");
        Ok(cmd)
    }

    pub fn exec(
        &self,
        program: impl AsRef<OsStr>,
        args: &[String],
        out: Option<&str>,
    ) -> Result<()> {
        let program = program.as_ref();
        eprintln!(
            "$ {} {}{}",
            program.to_string_lossy(),
            display_args(args),
            out.map(|p| format!(" > {p:?}")).unwrap_or_default()
        );
        if self.dry_run {
            return Ok(());
        }
        let mut cmd = self.command(program, args)?;
        // Write redirected output beside its destination and replace only on success.
        let pending = out.map(|p| {
            self.root
                .join(p)
                .with_extension(format!("tmp-{}", std::process::id()))
        });
        if let Some(path) = &pending {
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent)?;
            }
            cmd.stdout(OpenOptions::new().write(true).create_new(true).open(path)?);
        }
        let result = cmd
            .status()
            .with_context(|| format!("cannot start {}", program.to_string_lossy()));
        match result {
            Ok(status) if status.success() => {
                if let (Some(path), Some(dest)) = (pending, out) {
                    fs::rename(path, self.root.join(dest))?;
                }
                Ok(())
            }
            other => {
                if let Some(path) = pending {
                    let _ = fs::remove_file(path);
                }
                let status = other?;
                bail!("{} failed ({status})", program.to_string_lossy())
            }
        }
    }

    pub fn capture(&self, program: impl AsRef<OsStr>, args: &[String]) -> Result<String> {
        eprintln!(
            "$ {} {}",
            program.as_ref().to_string_lossy(),
            display_args(args)
        );
        if self.dry_run {
            return Ok(String::new());
        }
        let output = self
            .command(program.as_ref(), args)?
            .stderr(Stdio::inherit())
            .output()
            .with_context(|| format!("cannot start {}", program.as_ref().to_string_lossy()))?;
        if !output.status.success() {
            bail!(
                "{} failed ({})",
                program.as_ref().to_string_lossy(),
                output.status
            );
        }
        Ok(String::from_utf8(output.stdout)?)
    }

    pub fn cargo(&self, args: &[&str]) -> Result<()> {
        self.exec(cargo(), &strings(args), None)
    }

    pub fn run_args(&self, package: &str, bin: Option<&str>, args: &[&str]) -> Vec<String> {
        let mut result = strings(&["run", "--quiet", "--locked", "-p", package]);
        if self.release {
            result.push("--release".into());
        }
        if let Some(bin) = bin {
            result.extend(strings(&["--bin", bin]));
        }
        result.push("--".into());
        result.extend(strings(args));
        result
    }

    pub fn curve(&self, args: &[&str], out: Option<&str>) -> Result<()> {
        self.exec(cargo(), &self.run_args("curve_fitter", None, args), out)
    }

    pub fn curve_capture(&self, args: &[&str]) -> Result<String> {
        self.capture(cargo(), &self.run_args("curve_fitter", None, args))
    }

    pub fn python_path(&self, name: &str) -> Result<PathBuf> {
        let stem = name.trim_end_matches(".py").replace('-', "_");
        if stem.contains('/') || stem.contains('\\') || stem.starts_with('.') {
            bail!("use a script name from './run.sh python list'");
        }
        for prefix in ["scripts", "scripts/archive"] {
            let path = PathBuf::from(prefix).join(format!("{stem}.py"));
            if self.root.join(&path).is_file() {
                return Ok(path);
            }
        }
        bail!("unknown Python script '{name}'; use './run.sh python list'")
    }

    pub fn python(&self, name: &str, args: &[&str], out: Option<&str>) -> Result<()> {
        let path = self.python_path(name)?;
        let mut argv = vec![path.to_string_lossy().into_owned()];
        argv.extend(strings(args));
        self.exec(python(), &argv, out)
    }

    pub fn python_capture(&self, name: &str, args: &[&str]) -> Result<String> {
        let mut argv = vec![self.python_path(name)?.to_string_lossy().into_owned()];
        argv.extend(strings(args));
        self.capture(python(), &argv)
    }

    pub fn dirs(&self) -> Result<()> {
        if !self.dry_run {
            for dir in ["curves", "plots", "reports", "images", "data/milky_way"] {
                fs::create_dir_all(self.root.join(dir))?;
            }
        }
        Ok(())
    }

    pub fn env_default(&mut self, key: &str, value: &str) {
        if !self.env.contains_key(key) {
            self.env.insert(
                key.into(),
                std::env::var(key).unwrap_or_else(|_| value.into()),
            );
        }
    }

    pub fn cuda_build(&mut self) -> Result<()> {
        if self.cuda_ready {
            return Ok(());
        }
        let nvcc = std::env::var_os("HORB_NVCC").unwrap_or_else(|| {
            if Path::new("/usr/local/cuda/bin/nvcc").is_file() {
                "/usr/local/cuda/bin/nvcc".into()
            } else {
                "nvcc".into()
            }
        });
        let arch = std::env::var("HORB_CUDA_ARCH").unwrap_or_else(|_| "sm_87".into());
        let dest = self.root.join("cuda_kernels/cuda/libhorb_cuda.so");
        let pending = dest.with_extension(format!("so.tmp-{}", std::process::id()));
        let args = vec![
            format!("-arch={arch}"),
            "-O3".into(),
            "--compiler-options".into(),
            "-fPIC".into(),
            "-shared".into(),
            "cuda_kernels/cuda/horb_orbital_engine.cu".into(),
            "-o".into(),
            pending.to_string_lossy().into_owned(),
        ];
        let result = self.exec(nvcc, &args, None).context("CUDA build failed; install the CUDA toolkit or set HORB_NVCC; use 'build cpu' on hosts without CUDA");
        if !self.dry_run {
            if result.is_ok() {
                fs::rename(&pending, &dest)?;
            } else {
                let _ = fs::remove_file(&pending);
            }
        }
        result?;
        self.cuda_ready = true;
        Ok(())
    }

    pub fn cuda(&mut self, bin: &str, args: &[&str], out: Option<&str>) -> Result<()> {
        // Validate before attempting a potentially expensive build.
        if bin.contains('/')
            || !self
                .root
                .join(format!("cuda_kernels/src/bin/{bin}.rs"))
                .is_file()
        {
            bail!("unknown CUDA binary '{bin}'; use './run.sh run list'");
        }
        self.dirs()?;
        self.cuda_build()?;
        self.exec(
            cargo(),
            &self.run_args("cuda_kernels", Some(bin), args),
            out,
        )
    }

    pub fn write(&self, path: &str, text: &str) -> Result<()> {
        if !self.dry_run {
            fs::write(self.root.join(path), text).with_context(|| format!("write {path}"))?;
        }
        Ok(())
    }

    pub fn copy(&self, from: &str, to: &str) -> Result<()> {
        eprintln!("copy {from} -> {to}");
        if !self.dry_run {
            fs::copy(self.root.join(from), self.root.join(to))?;
        }
        Ok(())
    }
}

pub fn cargo() -> OsString {
    std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into())
}
pub fn python() -> OsString {
    std::env::var_os("HORB_PYTHON").unwrap_or_else(|| "python3".into())
}
fn display_args(args: &[String]) -> String {
    args.iter()
        .map(|s| {
            if s.chars()
                .all(|c| c.is_ascii_alphanumeric() || "-_/.:=,+".contains(c))
            {
                s.clone()
            } else {
                format!("{s:?}")
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;

    #[test]
    fn failed_redirect_preserves_previous_result_and_cleans_temp_file() {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("curve.csv"), "previous good result\n").unwrap();
        let r = Runner::new(dir.path().into(), false, false);
        assert!(r
            .exec(
                "sh",
                &strings(&["-c", "printf 'partial result'; exit 23"]),
                Some("curve.csv")
            )
            .is_err());
        assert_eq!(
            fs::read_to_string(dir.path().join("curve.csv")).unwrap(),
            "previous good result\n"
        );
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
        r.exec(
            "sh",
            &strings(&["-c", "printf 'complete result\\n'"]),
            Some("curve.csv"),
        )
        .unwrap();
        assert_eq!(
            fs::read_to_string(dir.path().join("curve.csv")).unwrap(),
            "complete result\n"
        );
    }

    #[test]
    fn arguments_are_literal_and_children_run_from_the_repository_root() {
        let dir = tempfile::tempdir().unwrap();
        let r = Runner::new(dir.path().into(), false, false);
        let value = "spaces ; $(touch INJECTED) `touch INJECTED`";
        assert_eq!(
            r.capture("printf", &strings(&["%s", value])).unwrap(),
            value
        );
        assert!(!dir.path().join("INJECTED").exists());
        assert_eq!(
            r.capture("pwd", &[]).unwrap().trim(),
            dir.path().to_str().unwrap()
        );
        let ld = r
            .capture("sh", &strings(&["-c", "printf '%s' \"$LD_LIBRARY_PATH\""]))
            .unwrap();
        assert!(ld.starts_with(dir.path().join("cuda_kernels/cuda").to_str().unwrap()));
    }
}
