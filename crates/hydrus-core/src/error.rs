use thiserror::Error;

#[derive(Error, Debug)]
pub enum HydrusError {
    #[error("Invalid project: {0}")]
    Invalid(String),
    #[error("Not supported: {0}")]
    Unsupported(String),
    #[error("Numerical error: {0}")]
    Numerical(String),
    #[error("I/O error: {0}")]
    Io(String),
    #[error("Parse error: {0}")]
    Parse(String),
}
