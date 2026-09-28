//! The wrapper that carries an [`Error`] across the FFI boundary as a built-in
//! Python exception.
//!
//! `PyErr` and [`Error`] are both foreign here, so the orphan rule blocks
//! `impl From<Error> for PyErr`. [`PyRisearchError`] owns that conversion instead:
//! `?` reaches it from either error type, and PyO3 accepts it as a return type
//! because it converts back into `PyErr`.

use pyo3::exceptions::{
    PyFileNotFoundError, PyOSError, PyPermissionError, PyRuntimeError, PyValueError,
};
use pyo3::PyErr;
use risearch::Error;

/// A `PyErr` in transit, convertible from either error type `?` can encounter.
pub struct PyRisearchError(PyErr);

impl From<Error> for PyRisearchError {
    fn from(err: Error) -> Self {
        let msg = err.to_string();
        Self(match err {
            Error::Io(io) => match io.kind() {
                std::io::ErrorKind::NotFound => PyFileNotFoundError::new_err(msg),
                std::io::ErrorKind::PermissionDenied => PyPermissionError::new_err(msg),
                _ => PyOSError::new_err(msg),
            },
            Error::Index(_) | Error::Dsm(_) | Error::Input(_) | Error::Config(_) => {
                PyValueError::new_err(msg)
            }
            Error::Output(_) => PyOSError::new_err(msg),
            // `Error` is #[non_exhaustive]: a variant added upstream lands here.
            _ => PyRuntimeError::new_err(msg),
        })
    }
}

impl From<PyErr> for PyRisearchError {
    fn from(err: PyErr) -> Self {
        Self(err)
    }
}

impl From<PyRisearchError> for PyErr {
    fn from(err: PyRisearchError) -> Self {
        err.0
    }
}

/// [`Result`](std::result::Result) for functions exposed to Python.
pub type Result<T> = std::result::Result<T, PyRisearchError>;
