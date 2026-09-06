//! Micro-benchmark for the `--thread 1` pool cost seen by in-process callers.
//!
//! Spawns `THREADS` worker threads, each running `CALLS_PER_THREAD`
//! `mafft_rs::run_from` invocations with `--thread 1` on a small DNA
//! fixture, and reports wall time plus per-call mean. Every output is
//! hashed and compared so a pool change that altered results would be
//! caught here too.
//!
//! Usage: `cargo run --release -p mafft-rs --example thread1_pool_bench [fixture] [threads] [calls_per_thread]`

use std::collections::hash_map::DefaultHasher;
use std::ffi::OsString;
use std::hash::{Hash, Hasher};
use std::path::PathBuf;
use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let fixture: PathBuf = args
        .next()
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/input_handling/dna_mixed_case.fa")
        });
    let threads: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(8);
    let calls: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(25);

    let argv: Vec<OsString> = [
        "mafft-rs", "--quiet", "--thread", "1", "--retree", "1", "--maxiterate", "0",
    ]
    .iter()
    .map(OsString::from)
    .chain(std::iter::once(fixture.clone().into_os_string()))
    .collect();

    // Warm-up / reference output.
    let mut reference = Vec::new();
    mafft_rs::run_from(argv.clone(), &mut reference).expect("reference run");
    let reference_hash = hash(&reference);

    let start = Instant::now();
    let handles: Vec<_> = (0..threads)
        .map(|_| {
            let argv = argv.clone();
            std::thread::spawn(move || {
                let mut hashes = Vec::with_capacity(calls);
                for _ in 0..calls {
                    let mut out = Vec::new();
                    mafft_rs::run_from(argv.clone(), &mut out).expect("run");
                    hashes.push(hash(&out));
                }
                hashes
            })
        })
        .collect();
    let mut total = 0usize;
    let mut mismatches = 0usize;
    for h in handles {
        for hv in h.join().expect("worker panicked") {
            total += 1;
            if hv != reference_hash {
                mismatches += 1;
            }
        }
    }
    let elapsed = start.elapsed();
    println!(
        "fixture={} threads={} calls={} wall={:.3}s per_call={:.2}ms mismatches={}",
        fixture.display(),
        threads,
        total,
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1000.0 / total as f64,
        mismatches
    );
    if mismatches != 0 {
        std::process::exit(2);
    }
}

fn hash(bytes: &[u8]) -> u64 {
    let mut h = DefaultHasher::new();
    bytes.hash(&mut h);
    h.finish()
}
