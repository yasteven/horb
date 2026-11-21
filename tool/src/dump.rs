//! Export the current working tree, including modified and untracked source.
use anyhow::{bail, Context, Result};
use std::collections::BTreeSet;
use std::fs::{self, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Component, Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

pub fn dump(root: &Path, output: Option<&Path>, data: bool, dry_run: bool) -> Result<PathBuf> {
    let now = SystemTime::now().duration_since(UNIX_EPOCH)?;
    let output = output.map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from(format!(
            "dumps/cat_all_horb_{}_{:09}.txt",
            now.as_secs(),
            now.subsec_nanos()
        ))
    });
    let output = if output.is_absolute() {
        output
    } else {
        root.join(output)
    };
    if output.exists() {
        bail!(
            "dump output already exists: {}; choose another filename",
            output.display()
        );
    }
    let paths = source_paths(root, data)?;
    if paths.is_empty() {
        bail!("no source files found under {}", root.display());
    }
    if dry_run {
        eprintln!("dump {} source files -> {}", paths.len(), output.display());
        return Ok(output);
    }
    // Read before creating the output: no accidental self-inclusion or partial dump
    // when a tracked source file is unreadable. Deleted files remain in git status.
    let mut contents = Vec::new();
    let mut skipped = Vec::new();
    for relative in paths {
        let bytes = fs::read(root.join(&relative))
            .with_context(|| format!("read {}", relative.display()))?;
        match String::from_utf8(bytes) {
            Ok(text) if !text.contains('\0') => contents.push((relative, text)),
            _ => skipped.push(relative),
        }
    }
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut file = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&output)?,
    );
    writeln!(
        file,
        "# HORB current source dump\n# generated_unix_seconds: {}\n# root: {}",
        now.as_secs(),
        root.display()
    )?;
    writeln!(
        file,
        "# branch: {}\n# commit: {}",
        git_text(root, &["branch", "--show-current"]),
        git_text(root, &["rev-parse", "HEAD"])
    )?;
    writeln!(
        file,
        "# Working-tree status (includes local changes):\n{}",
        git_text(root, &["status", "--short", "--untracked-files=all"])
    )?;
    writeln!(file, "# Contents are working-tree bytes, not git HEAD.\n# Files: {}\n# Data CSVs: {}\n# Binary/non-UTF8 files omitted: {}\n", contents.len(), data, skipped.len())?;
    for (path, text) in &contents {
        writeln!(file, "# {:>8} bytes  {}", text.len(), path.display())?;
    }
    for path in &skipped {
        writeln!(file, "# skipped binary/non-UTF8: {}", path.display())?;
    }
    for (path, text) in &contents {
        writeln!(
            file,
            "\n===== BEGIN FILE: {} =====\n# bytes: {} | lines: {}",
            path.display(),
            text.len(),
            text.lines().count()
        )?;
        file.write_all(text.as_bytes())?;
        writeln!(file, "\n===== END FILE: {} =====", path.display())?;
    }
    file.flush()?;
    eprintln!(
        "Wrote {} files ({} bytes): {}",
        contents.len(),
        fs::metadata(&output)?.len(),
        output.display()
    );
    Ok(output)
}

fn git_text(root: &Path, args: &[&str]) -> String {
    Command::new("git")
        .args(args)
        .current_dir(root)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unavailable (source export without Git metadata)".into())
}

fn source_paths(root: &Path, data: bool) -> Result<BTreeSet<PathBuf>> {
    let mut paths = BTreeSet::new();
    match Command::new("git")
        .args([
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "-z",
        ])
        .current_dir(root)
        .output()
    {
        Ok(output) if output.status.success() => {
            for raw in output.stdout.split(|b| *b == 0).filter(|p| !p.is_empty()) {
                let path =
                    PathBuf::from(std::str::from_utf8(raw).context("source path is not UTF-8")?);
                if eligible(&path, data) && regular_without_symlinks(root, &path) {
                    paths.insert(path);
                }
            }
        }
        _ => walk(root, Path::new(""), data, &mut paths)?,
    }
    Ok(paths)
}

fn walk(root: &Path, dir: &Path, data: bool, paths: &mut BTreeSet<PathBuf>) -> Result<()> {
    for entry in fs::read_dir(root.join(dir))? {
        let entry = entry?;
        let relative = dir.join(entry.file_name());
        let kind = entry.file_type()?;
        if excluded(&relative) || kind.is_symlink() {
            continue;
        }
        if kind.is_dir() {
            walk(root, &relative, data, paths)?;
        } else if kind.is_file() && eligible(&relative, data) {
            paths.insert(relative);
        }
    }
    Ok(())
}

fn regular_without_symlinks(root: &Path, path: &Path) -> bool {
    let mut current = root.to_path_buf();
    for part in path.components() {
        if !matches!(part, Component::Normal(_)) {
            return false;
        }
        current.push(part);
        match fs::symlink_metadata(&current) {
            Ok(m) if !m.file_type().is_symlink() => {}
            _ => return false,
        }
    }
    current.is_file()
}

fn excluded(path: &Path) -> bool {
    path.components().any(|c| {
        let name = c.as_os_str().to_string_lossy();
        matches!(
            name.as_ref(),
            ".git"
                | "target"
                | "dumps"
                | "curves"
                | "plots"
                | "reports"
                | "images"
                | ".venv"
                | "venv"
                | "node_modules"
                | "__pycache__"
                | ".ssh"
        ) || name == ".env"
            || name.starts_with(".env.")
    })
}

fn eligible(path: &Path, data: bool) -> bool {
    if excluded(path) {
        return false;
    }
    let name = path.file_name().and_then(|x| x.to_str()).unwrap_or("");
    if name.starts_with("cat_all_")
        || name.starts_with("patch_")
        || name.ends_with(".pem")
        || name.ends_with(".key")
    {
        return false;
    }
    if name.starts_with("LICENSE")
        || matches!(
            name,
            ".gitignore" | ".gitattributes" | "Makefile" | "Dockerfile"
        )
    {
        return true;
    }
    let ext = path.extension().and_then(|x| x.to_str()).unwrap_or("");
    matches!(
        ext,
        "rs" | "cu"
            | "cuh"
            | "c"
            | "h"
            | "cpp"
            | "hpp"
            | "py"
            | "sh"
            | "toml"
            | "lock"
            | "md"
            | "yml"
            | "yaml"
            | "json"
    ) || (data && path.starts_with("data") && matches!(ext, "csv" | "dat" | "mrt" | "txt"))
        || name.starts_with("requirements") && ext == "txt"
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dump_reads_modified_and_untracked_files_and_refuses_overwrite() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path();
        assert!(Command::new("git")
            .args(["init", "--quiet"])
            .current_dir(root)
            .status()
            .unwrap()
            .success());
        fs::write(root.join("tracked.rs"), "old source\n").unwrap();
        assert!(Command::new("git")
            .args(["add", "tracked.rs"])
            .current_dir(root)
            .status()
            .unwrap()
            .success());
        fs::write(root.join("tracked.rs"), "current local changes\n").unwrap();
        fs::write(root.join("new.rs"), "new untracked code\n").unwrap();
        fs::write(root.join(".env"), "DO_NOT_EXPORT=secret\n").unwrap();
        fs::create_dir(root.join("data")).unwrap();
        fs::write(root.join("data/curve.csv"), "r,v\n1,2\n").unwrap();
        let out = dump(root, None, false, false).unwrap();
        let text = fs::read_to_string(&out).unwrap();
        assert!(text.contains("current local changes"));
        assert!(text.contains("new untracked code"));
        assert!(!text.contains("old source"));
        assert!(!text.contains("DO_NOT_EXPORT"));
        assert!(!text.contains("BEGIN FILE: data/curve.csv"));
        assert!(dump(root, Some(&out), false, false).is_err());
        let second = dump(root, None, true, false).unwrap();
        assert_ne!(out, second);
        assert!(fs::read_to_string(second)
            .unwrap()
            .contains("BEGIN FILE: data/curve.csv"));
    }

    #[test]
    fn source_export_fallback_skips_symlinks_and_dry_run_writes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("main.rs"), "fn main() {}\n").unwrap();
        #[cfg(unix)]
        std::os::unix::fs::symlink("/etc/passwd", dir.path().join("outside.rs")).unwrap();
        let planned = dump(dir.path(), None, false, true).unwrap();
        assert!(!planned.exists());
        assert!(!dir.path().join("dumps").exists());
        let out = dump(dir.path(), None, false, false).unwrap();
        let text = fs::read_to_string(out).unwrap();
        assert!(text.contains("fn main() {}"));
        assert!(!text.contains("outside.rs"));
    }
    #[test]
    fn includes_source_but_excludes_generated_and_credentials() {
        for p in [
            "tool/src/main.rs",
            "cuda_kernels/cuda/kernel.cu",
            "Cargo.lock",
            ".gitignore",
            "README.md",
            "requirements.txt",
        ] {
            assert!(eligible(Path::new(p), false), "{p}");
        }
        for p in [
            "target/debug/x.rs",
            "dumps/old.rs",
            ".env",
            ".env.local",
            ".ssh/config.json",
            "plots/data.json",
            "cat_all_horb.txt",
            "key.pem",
        ] {
            assert!(!eligible(Path::new(p), true), "{p}");
        }
        assert!(!eligible(Path::new("data/milky_way/curve.csv"), false));
        assert!(eligible(Path::new("data/milky_way/curve.csv"), true));
    }
}
