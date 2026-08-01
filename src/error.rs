use std::fmt;

/// Errors that can occur during matrix operations.
///
/// Returned via `Result` when the input is structurally well-formed (e.g.
/// square or dimension-matched) but mathematically unsuitable for the
/// operation (e.g. singular for inversion, not positive definite for
/// Cholesky).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixError {
    /// The matrix is singular (non-invertible).
    ///
    /// Occurs in operations such as inversion, solving linear systems,
    /// and LU decomposition when the matrix has zero determinant.
    Singular,

    /// The matrix is not positive definite.
    ///
    /// Occurs in Cholesky (LDL^T) decomposition when the matrix is symmetric
    /// but has a non-positive pivot, meaning it is not positive definite.
    NotPositiveDefinite,
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixError::Singular => write!(f, "Error: The matrix is singular."),
            MatrixError::NotPositiveDefinite => {
                write!(f, "Error: The matrix is not positive definite.")
            }
        }
    }
}

impl std::error::Error for MatrixError {}
