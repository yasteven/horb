// data_io/examples/inspect_rc.rs

use anyhow::Result;
use data_io::RotationCurve;

fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("usage: cargo run -p data_io --example inspect_rc -- <csv>");

    let rc = RotationCurve::from_csv(path)?;

    println!("name: {}", rc.name);
    println!("rows: {}", rc.len());
    println!("r_min_kpc: {:.6}", rc.r_min());
    println!("r_max_kpc: {:.6}", rc.r_max());
    println!("has_baryons: {}", rc.has_baryons());

    Ok(())
}
