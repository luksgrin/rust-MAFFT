use std::path::PathBuf;

fn main() {
    let core_dir: PathBuf = ["../../mafft-upstream/core"].iter().collect();

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

    let mut source_paths: Vec<PathBuf> = library_sources
        .iter()
        .map(|f| core_dir.join(f))
        .collect();

    // Local wrappers for symbols stuck behind `static` in splittbfast.c
    // or living in C files with their own `main()` (addsingle.c,
    // disttbfast.c). Bodies are copies of the canonical implementations.
    let wrappers_dir: PathBuf = ["wrappers"].iter().collect();
    source_paths.push(wrappers_dir.join("parttree_helpers.c"));
    source_paths.push(wrappers_dir.join("msalignmm_instr.c"));
    source_paths.push(wrappers_dir.join("blend_helpers.c"));

    let mut build = cc::Build::new();
    build
        .files(&source_paths)
        .include(&core_dir)
        .include(&wrappers_dir)
        .opt_level(3)
        .std("c99")
        .define("enablemultithread", None)
        .warnings(false);

    // GCC 14 turned several long-standing warnings into errors that `-w`
    // does not silence. With GCC 14.4 / 15.2 / 15.3 (`-std=c99`) the
    // upstream sources hit exactly one of them:
    //   iteration.c:148, 338: passing argument 2 of 'writePre' from
    //   incompatible pointer type [-Wincompatible-pointer-types]
    // (`char (*)[256]` passed as `char **`; the upstream Makefile never
    // compiles iteration.c, so `make -C mafft-upstream/core` is unaffected).
    // The other GCC 14 / Clang 16 promotions are downgraded too so a
    // header or compiler change does not reopen this. Each flag is probed,
    // so compilers that do not know it (older GCC, MSVC) never see it.
    for flag in [
        "-Wno-error=incompatible-pointer-types",
        "-Wno-error=implicit-function-declaration",
        "-Wno-error=implicit-int",
        "-Wno-error=int-conversion",
        "-Wno-error=return-mismatch",
        "-Wno-error=incompatible-function-pointer-types",
    ] {
        build.flag_if_supported(flag);
    }

    // Platform-specific flags
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if target_os != "windows" {
        build.flag("-pthread");
    }

    build.compile("mafft_c");

    // GCC on macOS (e.g. Homebrew gcc-15) lowers `__thread`, used under
    // -Denablemultithread, to emulated TLS. rustc links with the system
    // clang, which does not provide `___emutls_get_address`, so pull in
    // GCC's own runtime for it.
    let compiler = build.get_compiler();
    if target_os == "macos" && compiler.is_like_gnu() && !compiler.is_like_clang() {
        let emutls_dir = std::process::Command::new(compiler.path())
            .arg("-print-file-name=libemutls_w.a")
            .output()
            .ok()
            .map(|out| PathBuf::from(String::from_utf8_lossy(&out.stdout).trim()))
            .filter(|lib| lib.is_absolute() && lib.exists())
            .and_then(|lib| lib.parent().map(PathBuf::from));
        if let Some(dir) = emutls_dir {
            println!("cargo:rustc-link-search=native={}", dir.display());
            println!("cargo:rustc-link-lib=static=emutls_w");
        }
    }

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
