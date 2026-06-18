use std::env;
use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let cuda_dir = manifest_dir.join("cuda");

    println!("cargo:rustc-link-search=native={}", cuda_dir.display());
    println!("cargo:rustc-link-lib=dylib=horb_cuda");
}
