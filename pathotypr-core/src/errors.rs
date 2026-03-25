//! src/errors.rs
//!
//! Module for custom error handling throughout the application.
//!
//! Defines an `AppError` enum that unifies the different types of errors that can
//! occur, such as I/O, parsing, model errors, etc. This allows for cleaner
//! and more specific error handling than using `anyhow`.

use indicatif::style::TemplateError;
use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

/// A type alias for the application's standard Result type.
pub type AppResult<T> = Result<T, AppError>;

/// A token that can be used to check if a task should be cancelled.
/// Clone is cheap as it only clones an Arc.
#[derive(Clone, Default)]
pub struct CancellationToken {
    cancelled: Arc<AtomicBool>,
}

impl std::fmt::Debug for CancellationToken {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CancellationToken")
            .field("cancelled", &self.cancelled.load(Ordering::SeqCst))
            .finish()
    }
}

impl CancellationToken {
    /// Create a new cancellation token.
    pub fn new() -> Self {
        Self {
            cancelled: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Request cancellation.
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::SeqCst);
    }

    /// Check if cancellation has been requested.
    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::SeqCst)
    }

    /// Reset the cancellation state.
    pub fn reset(&self) {
        self.cancelled.store(false, Ordering::SeqCst);
    }

    /// Check if cancelled and return an error if so.
    pub fn check(&self) -> AppResult<()> {
        if self.is_cancelled() {
            Err(AppError::Cancelled)
        } else {
            Ok(())
        }
    }
}

/// Helper function to check an optional cancellation token.
/// Returns Ok(()) if token is None or not cancelled, Err(Cancelled) if cancelled.
#[inline]
pub fn check_cancelled(token: &Option<CancellationToken>) -> AppResult<()> {
    if let Some(t) = token {
        t.check()
    } else {
        Ok(())
    }
}

/// Helper struct for managing cancellation in parallel processing.
/// Combines an AtomicBool flag (for fast checking across threads) with
/// an optional CancellationToken (for external cancellation requests).
pub struct ParallelCancellation<'a> {
    /// Fast flag that can be checked without locking.
    cancelled: AtomicBool,
    /// External cancellation token.
    token: &'a Option<CancellationToken>,
}

impl<'a> ParallelCancellation<'a> {
    /// Create a new parallel cancellation helper.
    pub fn new(token: &'a Option<CancellationToken>) -> Self {
        Self {
            cancelled: AtomicBool::new(false),
            token,
        }
    }

    /// Check if cancellation has been requested.
    /// This is optimized for use in hot loops - first checks the fast flag,
    /// then checks the token if the flag isn't set.
    #[inline]
    pub fn is_cancelled(&self) -> bool {
        // Fast path: check the local flag first
        if self.cancelled.load(Ordering::Relaxed) {
            return true;
        }
        // Slow path: check the token and update the flag if cancelled
        if let Some(token) = self.token {
            if token.is_cancelled() {
                self.cancelled.store(true, Ordering::Relaxed);
                return true;
            }
        }
        false
    }

    /// Check if cancelled and return Err if so. For use after parallel processing.
    pub fn check_after(&self) -> AppResult<()> {
        check_cancelled(self.token)
    }
}

/// An enum representing all possible errors in the application.
#[derive(Debug)]
pub enum AppError {
    /// Error related to input/output operations.
    Io(std::io::Error),
    /// Error when deserializing or serializing data (e.g., bincode).
    Serialization(Box<bincode::ErrorKind>),
    /// Error when parsing data (e.g., CSV, FASTA).
    Parsing(String),
    /// Error when there is not enough data to train or predict.
    NotEnoughData(String),
    /// Error from the progress bar template.
    Template(TemplateError),
    /// Error when task is cancelled by user.
    Cancelled,
    /// A generic error for other cases.
    Generic(String),
}

// Display implementation to show user-friendly errors.
impl fmt::Display for AppError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AppError::Io(err) => write!(f, "IO Error: {}", err),
            AppError::Serialization(err) => write!(f, "Serialization Error: {}", err),
            AppError::Parsing(msg) => write!(f, "Parsing Error: {}", msg),
            AppError::NotEnoughData(msg) => write!(f, "Not Enough Data: {}", msg),
            AppError::Template(err) => write!(f, "UI Template Error: {}", err),
            AppError::Cancelled => write!(f, "Task cancelled by user"),
            AppError::Generic(msg) => write!(f, "Error: {}", msg),
        }
    }
}

// Implementation of std::error::Error.
impl std::error::Error for AppError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            AppError::Io(err) => Some(err),
            AppError::Serialization(err) => Some(err),
            AppError::Template(err) => Some(err),
            AppError::Cancelled => None,
            _ => None,
        }
    }
}

// Conversions from other error types to AppError.

impl From<std::io::Error> for AppError {
    fn from(err: std::io::Error) -> Self {
        AppError::Io(err)
    }
}

impl From<Box<bincode::ErrorKind>> for AppError {
    fn from(err: Box<bincode::ErrorKind>) -> Self {
        AppError::Serialization(err)
    }
}

impl From<csv::Error> for AppError {
    fn from(err: csv::Error) -> Self {
        AppError::Parsing(format!("CSV error: {}", err))
    }
}

impl From<TemplateError> for AppError {
    fn from(err: TemplateError) -> Self {
        AppError::Template(err)
    }
}
