extern crate cbindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let package_name = env::var("CARGO_PKG_NAME").unwrap();
    
    // Path for the header file relative to the CARGO_MANIFEST_DIR
    // This will be <crate_root>/target/include/<package_name>.h
    let output_file = PathBuf::from(&crate_dir) // Start from crate_dir
        .join("target")
        .join("include")
        .join(format!("{}.h", package_name));

    // Create the target/include directory if it doesn't exist
    // This needs to be done relative to the crate_dir as well
    if let Some(parent_dir) = output_file.parent() {
        if !parent_dir.exists() {
            std::fs::create_dir_all(parent_dir)
                .expect("Failed to create target/include directory");
        }
    }


    match cbindgen::generate(&crate_dir) {
        Ok(bindings) => {
            bindings.write_to_file(&output_file);
        }
        Err(err) => {
            // Print the error to stderr for better diagnostics in cargo output
            eprintln!("Unable to generate bindings: {:?}", err);
            panic!("Unable to generate bindings: {:?}", err);
        }
    }
}
