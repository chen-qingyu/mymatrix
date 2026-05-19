use std::fmt;

/// Errors that can occur during matrix operations.
///
/// Returned via `Result` for scenarios where the input is mathematically
/// valid but the operation cannot produce a meaningful result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixError {
    /// The matrix is singular (non-invertible).
    ///
    /// Occurs in operations such as inversion, solving linear systems,
    /// and LU decomposition when the matrix has zero determinant.
    Singular,
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixError::Singular => write!(f, "Error: The matrix is singular."),
        }
    }
}

impl std::error::Error for MatrixError {}
