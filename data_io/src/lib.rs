// data_io/src/lib.rs

use anyhow::{Context, Result};
use serde::Deserialize;
use std::path::Path;

#[derive(Debug, Clone, Deserialize)]
pub struct RotationCurveRow {
    #[serde(rename = "R_kpc")]
    pub r_kpc: f64,

    #[serde(rename = "Vobs_kms")]
    pub vobs_kms: f64,

    #[serde(rename = "eV_kms")]
    pub ev_kms: f64,

    #[serde(rename = "Vgas_kms", default)]
    pub vgas_kms: Option<f64>,

    #[serde(rename = "Vdisk_kms", default)]
    pub vdisk_kms: Option<f64>,

    #[serde(rename = "Vbul_kms", default)]
    pub vbul_kms: Option<f64>,

    #[serde(default)]
    pub source: Option<String>,
}

#[derive(Debug, Clone)]
pub struct RotationCurve {
    pub name: String,
    pub rows: Vec<RotationCurveRow>,
}

impl RotationCurve {
    pub fn from_csv(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();

        let mut rdr = csv::Reader::from_path(path)
            .with_context(|| format!("failed to open rotation curve CSV: {}", path.display()))?;

        let mut rows = Vec::new();

        for rec in rdr.deserialize() {
            let row: RotationCurveRow =
                rec.with_context(|| format!("failed to parse row in {}", path.display()))?;
            rows.push(row);
        }

        if rows.is_empty() {
            anyhow::bail!("rotation curve CSV has no rows: {}", path.display());
        }

        rows.sort_by(|a, b| a.r_kpc.partial_cmp(&b.r_kpc).unwrap());

        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("rotation_curve")
            .to_string();

        Ok(Self { name, rows })
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn r_min(&self) -> f64 {
        self.rows.first().map(|r| r.r_kpc).unwrap_or(0.0)
    }

    pub fn r_max(&self) -> f64 {
        self.rows.last().map(|r| r.r_kpc).unwrap_or(0.0)
    }

    pub fn has_baryons(&self) -> bool {
        self.rows.iter().any(|r| {
            r.vgas_kms.unwrap_or(0.0) != 0.0
                || r.vdisk_kms.unwrap_or(0.0) != 0.0
                || r.vbul_kms.unwrap_or(0.0) != 0.0
        })
    }
}
