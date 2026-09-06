//! Python bindings for MAFFT-rs multiple sequence alignment.
//!
//! This crate is the native half of the `pymafft` package. It owns the
//! result types (`AlignmentResult`, `AlignedSequence`), the `MafftError`
//! exception, the input duck-typing, and three private entry points
//! (`_run`, `_run_file`, `_run_fasta_string`) that hand a CLI flag list to
//! [`mafft_rs::run_from_seqs`]. The public functions `align`, `align_file`,
//! `align_fasta_string` and `run` live in `python/pymafft/__init__.py`, where
//! every keyword option is turned into the corresponding `mafft-rs` flag.
//!
//! Routing through `run_from_seqs` means the Python package never
//! re-derives what a flag means: `--auto`'s size heuristic,
//! `--adjustdirection`'s strand detection, `--nuc` / `--amino`'s type
//! forcing and case fold, `--thread`'s pool, `--reorder`, the scoring
//! matrices and the gap penalties are all decided by the same code the
//! command line runs, so the Python and CLI outputs are byte-identical.
//!
//! ```python
//! import pymafft
//!
//! result = pymafft.align(["ACDEFGHIK", "ACDEFHIK", "ACDHIK"])
//! for seq in result:
//!     print(f"{seq.name}: {seq.sequence}")
//!
//! # The panaroo flag set, over in-memory records:
//! result = pymafft.align(records, strategy="auto", adjust_direction=True,
//!                        threads=1, seq_type="nuc")
//! ```

use std::sync::Mutex;

use pyo3::create_exception;
use pyo3::exceptions::{PyImportError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyList};

use mafft_io::{read_fasta, read_fasta_from_reader};
use mafft_rs::{run_from_seqs, MultipleAlignment, Progress, SeqType, Sequence, SequenceSet, StderrProgress};

create_exception!(
    pymafft,
    MafftError,
    PyValueError,
    "Raised when the alignment fails.\n\n\
     Carries the exact stderr text and exit code the `mafft-rs` command line\n\
     would have produced for the same input and flags: `.message` is the\n\
     text (also `str(err)`), `.code` the exit code. A `ValueError` subclass,\n\
     so `except ValueError` keeps working."
);

/// A single aligned sequence with name and gapped data.
#[pyclass]
#[derive(Clone)]
struct AlignedSequence {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    sequence: String,
}

#[pymethods]
impl AlignedSequence {
    fn __repr__(&self) -> String {
        format!("AlignedSequence(name='{}', len={})", self.name, self.sequence.len())
    }

    fn __str__(&self) -> String {
        format!(">{}\n{}", self.name, self.sequence)
    }

    /// Return the sequence without gap characters.
    fn ungapped(&self) -> String {
        self.sequence.chars().filter(|&c| c != '-').collect()
    }

    /// Length of the aligned sequence (including gaps).
    fn __len__(&self) -> usize {
        self.sequence.len()
    }
}

/// Result of a multiple sequence alignment.
#[pyclass]
#[derive(Clone)]
struct AlignmentResult {
    #[pyo3(get)]
    sequences: Vec<AlignedSequence>,
    #[pyo3(get)]
    width: usize,
    #[pyo3(get)]
    score: f64,
}

#[pymethods]
impl AlignmentResult {
    fn __repr__(&self) -> String {
        format!(
            "AlignmentResult(nseq={}, width={}, score={:.1})",
            self.sequences.len(),
            self.width,
            self.score,
        )
    }

    fn __len__(&self) -> usize {
        self.sequences.len()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<AlignmentResultIter>> {
        let iter = AlignmentResultIter {
            inner: slf.sequences.clone(),
            index: 0,
        };
        Py::new(slf.py(), iter)
    }

    /// Get a sequence by index.
    fn __getitem__(&self, idx: usize) -> PyResult<AlignedSequence> {
        self.sequences.get(idx).cloned().ok_or_else(|| {
            PyValueError::new_err(format!("index {} out of range (nseq={})", idx, self.sequences.len()))
        })
    }

    /// Return the alignment as a FASTA-formatted string.
    fn to_fasta(&self) -> String {
        self.sequences
            .iter()
            .map(|s| format!(">{}\n{}", s.name, s.sequence))
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// Return the alignment as a list of (name, sequence) tuples.
    fn to_tuples(&self) -> Vec<(String, String)> {
        self.sequences
            .iter()
            .map(|s| (s.name.clone(), s.sequence.clone()))
            .collect()
    }

    /// Return the alignment as a `Bio.Align.MultipleSeqAlignment`.
    ///
    /// Each row becomes a `Bio.SeqRecord.SeqRecord` whose `id` is the
    /// sequence name and whose `description` is empty. Biopython is
    /// imported lazily, here and only here; it is not a dependency of
    /// pymafft. Raises `ImportError` if it is not installed.
    fn to_biopython<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let import = |module: &str, name: &str| -> PyResult<Bound<'py, PyAny>> {
            py.import(module)
                .map_err(|e| {
                    PyImportError::new_err(format!(
                        "AlignmentResult.to_biopython() requires Biopython \
                         (`pip install biopython`): {e}"
                    ))
                })?
                .getattr(name)
        };
        let seq_cls = import("Bio.Seq", "Seq")?;
        let record_cls = import("Bio.SeqRecord", "SeqRecord")?;
        let msa_cls = import("Bio.Align", "MultipleSeqAlignment")?;

        let records = PyList::empty(py);
        for s in &self.sequences {
            let seq = seq_cls.call1((s.sequence.as_str(),))?;
            let kwargs = [("id", s.name.as_str()), ("description", "")].into_py_dict(py)?;
            records.append(record_cls.call((seq,), Some(&kwargs))?)?;
        }
        msa_cls.call1((records,))
    }

    /// Number of sequences.
    #[getter]
    fn nseq(&self) -> usize {
        self.sequences.len()
    }
}

#[pyclass]
struct AlignmentResultIter {
    inner: Vec<AlignedSequence>,
    index: usize,
}

#[pymethods]
impl AlignmentResultIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<AlignedSequence> {
        if self.index < self.inner.len() {
            let item = self.inner[self.index].clone();
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}

// ---------------------------------------------------------------------------
// Input conversion
// ---------------------------------------------------------------------------

/// Try to convert a Python object with `.id` and `.seq` attributes (duck-
/// types `Bio.SeqRecord.SeqRecord` and similar) into a `(name, seq)` pair.
/// Returns `None` if either attribute is missing.
fn try_record_to_pair(obj: &Bound<'_, PyAny>) -> Option<(String, String)> {
    let id_obj = obj.getattr("id").ok()?;
    let seq_obj = obj.getattr("seq").ok()?;
    let id: String = id_obj.extract().ok()?;
    // `seq` may be a `Bio.Seq.Seq` (custom type) — `str()` always works.
    let seq_str = seq_obj.str().ok()?;
    let seq: String = seq_str.extract().ok()?;
    Some((id, seq))
}

/// Accept, in order of try: a list of strings (auto-named `seq_N`), a list
/// of `(name, seq)` tuples, or an iterable of SeqRecord-like objects.
fn extract_pairs(sequences: &Bound<'_, PyAny>) -> PyResult<Vec<(String, String)>> {
    if let Ok(str_list) = sequences.extract::<Vec<String>>() {
        return Ok(str_list
            .into_iter()
            .enumerate()
            .map(|(i, s)| (format!("seq_{}", i + 1), s))
            .collect());
    }
    if let Ok(tuple_list) = sequences.extract::<Vec<(String, String)>>() {
        return Ok(tuple_list);
    }
    if let Ok(iter) = sequences.try_iter() {
        let mut pairs: Vec<(String, String)> = Vec::new();
        for (i, item) in iter.enumerate() {
            let item = item?;
            match try_record_to_pair(&item) {
                Some(p) => pairs.push(p),
                None => {
                    return Err(PyValueError::new_err(format!(
                        "sequences[{}] is not a string, (name, seq) tuple, or \
                         object with .id and .seq attributes",
                        i
                    )))
                }
            }
        }
        return Ok(pairs);
    }
    Err(PyValueError::new_err(
        "sequences must be a list of strings, (name, sequence) tuples, \
         or SeqRecord-like objects",
    ))
}

/// Build the `SequenceSet` handed to the flag layer. The type is left
/// `Unknown` so `run_from_seqs` detects it (or honours `--nuc` / `--amino`)
/// exactly as the command line does for a FASTA file, and the residues are
/// passed raw so its normalisation and case fold apply.
fn build_input(pairs: Vec<(String, String)>) -> PyResult<SequenceSet> {
    if pairs.is_empty() {
        return Err(PyValueError::new_err("no sequences provided"));
    }
    Ok(SequenceSet {
        sequences: pairs
            .into_iter()
            .map(|(name, data)| Sequence { name, data: data.into_bytes() })
            .collect(),
        seq_type: SeqType::Unknown,
    })
}

// ---------------------------------------------------------------------------
// Progress sink
// ---------------------------------------------------------------------------

/// Forwards each progress line to a Python callable.
///
/// The alignment runs with the GIL released, so every message re-acquires
/// it for the duration of the call. An exception raised by the callable is
/// stored and re-raised once the run returns; later messages are dropped.
struct PyProgress {
    callback: Py<PyAny>,
    error: Mutex<Option<PyErr>>,
}

impl Progress for PyProgress {
    fn message(&self, msg: &str) {
        let mut slot = match self.error.lock() {
            Ok(slot) => slot,
            Err(poisoned) => poisoned.into_inner(),
        };
        if slot.is_some() {
            return;
        }
        Python::with_gil(|py| {
            if let Err(e) = self.callback.call1(py, (msg,)) {
                *slot = Some(e);
            }
        });
    }
}

// ---------------------------------------------------------------------------
// Running
// ---------------------------------------------------------------------------

fn mafft_error_to_py(py: Python<'_>, err: mafft_rs::MafftError) -> PyErr {
    let py_err = MafftError::new_err(err.message().to_string());
    // Attach the CLI's exit code and message as attributes. `str(err)` is
    // already the message because it is the sole constructor argument.
    let value = py_err.value(py);
    let _ = value.setattr("code", err.code());
    let _ = value.setattr("message", err.message());
    py_err
}

fn msa_to_result(msa: MultipleAlignment) -> AlignmentResult {
    let width = msa.width();
    let sequences: Vec<AlignedSequence> = msa
        .sequences
        .into_iter()
        .zip(msa.names)
        .map(|(seq, name)| AlignedSequence {
            name,
            sequence: String::from_utf8_lossy(&seq).into_owned(),
        })
        .collect();
    AlignmentResult { width, score: msa.score, sequences }
}

/// Run `mafft-rs <args>` on `input` in-process, through the same flag layer
/// as the command line, with the GIL released for the duration.
fn run_with_flags(
    py: Python<'_>,
    args: Vec<String>,
    input: &SequenceSet,
    progress: Option<Py<PyAny>>,
) -> PyResult<AlignmentResult> {
    let mut argv = Vec::with_capacity(args.len() + 1);
    argv.push("mafft-rs".to_string());
    argv.extend(args);

    let outcome = match progress {
        Some(callback) => {
            let sink = PyProgress { callback, error: Mutex::new(None) };
            let outcome = py.allow_threads(|| run_from_seqs(argv, input, &sink));
            let callback_error = sink.error.into_inner().unwrap_or_else(|p| p.into_inner());
            if let Some(e) = callback_error {
                return Err(e);
            }
            outcome
        }
        // Without a sink, messages go to stderr exactly as on the command
        // line (`--quiet` in `args` silences them there too).
        None => py.allow_threads(|| run_from_seqs(argv, input, &StderrProgress)),
    };

    outcome.map(msa_to_result).map_err(|e| mafft_error_to_py(py, e))
}

/// Align in-memory sequences with an explicit `mafft-rs` flag list.
///
/// `args` are the command-line flags without the program name and without
/// an INPUT path (the sequences take its place). `sequences` accepts the
/// shapes documented for `pymafft.align`. `progress` is an optional
/// callable receiving each progress line; when `None`, progress goes to
/// stderr like the CLI (pass `--quiet` in `args` to silence it).
#[pyfunction]
#[pyo3(signature = (args, sequences, progress=None))]
fn _run(
    py: Python<'_>,
    args: Vec<String>,
    sequences: &Bound<'_, PyAny>,
    progress: Option<Py<PyAny>>,
) -> PyResult<AlignmentResult> {
    let input = build_input(extract_pairs(sequences)?)?;
    run_with_flags(py, args, &input, progress)
}

/// Read a FASTA file with the CLI's reader and align it with `args`.
#[pyfunction]
#[pyo3(signature = (args, path, progress=None))]
fn _run_file(
    py: Python<'_>,
    args: Vec<String>,
    path: &str,
    progress: Option<Py<PyAny>>,
) -> PyResult<AlignmentResult> {
    let input = read_fasta(path).map_err(|e| {
        PyValueError::new_err(format!("error reading FASTA file '{}': {}", path, e))
    })?;
    run_with_flags(py, args, &input, progress)
}

/// Parse a FASTA string with the CLI's reader and align it with `args`.
#[pyfunction]
#[pyo3(signature = (args, fasta_string, progress=None))]
fn _run_fasta_string(
    py: Python<'_>,
    args: Vec<String>,
    fasta_string: &str,
    progress: Option<Py<PyAny>>,
) -> PyResult<AlignmentResult> {
    let reader = std::io::Cursor::new(fasta_string.as_bytes());
    let input = read_fasta_from_reader(std::io::BufReader::new(reader)).map_err(|e| {
        PyValueError::new_err(format!("error parsing FASTA string: {}", e))
    })?;
    run_with_flags(py, args, &input, progress)
}

/// pymafft: Python bindings for MAFFT-rs multiple sequence alignment.
#[pymodule]
fn pymafft(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<AlignedSequence>()?;
    m.add_class::<AlignmentResult>()?;
    m.add("MafftError", m.py().get_type::<MafftError>())?;
    m.add_function(wrap_pyfunction!(_run, m)?)?;
    m.add_function(wrap_pyfunction!(_run_file, m)?)?;
    m.add_function(wrap_pyfunction!(_run_fasta_string, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
