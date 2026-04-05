use std::path::PathBuf;

fn main() {
    let core_dir: PathBuf = ["../../core"].iter().collect();

    // Library source files (no main function).
    // These are the shared objects used across MAFFT binaries.
    let library_sources = [
        "mtxutl.c",
        "io.c",
        "mltaln9.c",
        "tddis.c",
        "constants.c",
        "defs.c",
        "Galign11.c",
        "Lalign11.c",
        "genalign11.c",
        "fft.c",
        "fftFunctions.c",
        "Salignmm.c",
        "Dalignmm.c",
        "partSalignmm.c",
        "Lalignmm.c",
        "MSalignmm.c",
        "MSalign11.c",
        "Falign.c",
        "Falign_localhom.c",
        "SAalignmm.c",
        "rna.c",
        "nj.c",
        "treeOperation.c",
        "tditeration.c",
        "addfunctions.c",
        "pairlocalalign.c",
        "suboptalign11.c",
        "iteration.c",
    ];

    let source_paths: Vec<PathBuf> = library_sources
        .iter()
        .map(|f| core_dir.join(f))
        .collect();

    let mut build = cc::Build::new();
    build
        .files(&source_paths)
        .include(&core_dir)
        .opt_level(3)
        .std("c99")
        .define("enablemultithread", None)
        .warnings(false);

    // Platform-specific flags
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if target_os != "windows" {
        build.flag("-pthread");
    }

    build.compile("mafft_c");

    // Link system libraries
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=pthread");

    // Re-run if any C source changes
    for src in &source_paths {
        println!("cargo:rerun-if-changed={}", src.display());
    }
    for header in &["mltaln.h", "mtxutl.h", "mafft.h", "fft.h", "dp.h", "functions.h"] {
        println!("cargo:rerun-if-changed={}", core_dir.join(header).display());
    }
}
