//! src/errors.rs
//!
//! Module for custom error handling throughout the application.
//!
//! Defines an `AppError` enum that unifies the different types of errors that can
//! occur, such as I/O, parsing, model errors, etc. This allows for cleaner
//! and more specific error handling than using `anyhow`.

use indicatif::style::TemplateError;
use std::fmt;

/// A type alias for the application's standard Result type.
pub type AppResult<T> = Result<T, AppError>;

/// An enum representing all possible errors in the application.
#[derive(Debug)]
pub enum AppError {
    /// Error related to input/output operations.
    Io(std::io::Error),
    /// Error when deserializing or serializing data (e.g., bincode).
    Serialization(Box<bincode::ErrorKind>),
    /// Error when parsing data (e.g., CSV, FASTA).
    Parsing(String),
    /// Error originating from the machine learning library.
    SmartCore(smartcore::error::Failed),
    /// Error when there is not enough data to train or predict.
    NotEnoughData(String),
    /// Error from the progress bar template.
    Template(TemplateError),
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
            AppError::SmartCore(err) => write!(f, "Model Error: {}", err),
            AppError::NotEnoughData(msg) => write!(f, "Not Enough Data: {}", msg),
            AppError::Template(err) => write!(f, "UI Template Error: {}", err),
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
            AppError::SmartCore(err) => Some(err),
            AppError::Template(err) => Some(err),
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

impl From<smartcore::error::Failed> for AppError {
    fn from(err: smartcore::error::Failed) -> Self {
        AppError::SmartCore(err)
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
