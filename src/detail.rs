use crate::Matrix;

/// Check whether the index is valid (`begin <= pos < end`).
///
/// # Panics
/// Panics if `pos` is out of the range [`begin`, `end`).
#[inline]
pub fn check_bounds(pos: usize, begin: usize, end: usize) {
    if pos < begin || pos >= end {
        panic!("Error: Index out of range.");
    }
}

/// Check whether the container is not empty.
///
/// # Panics
/// Panics if `size` is 0.
#[inline]
pub fn check_empty(size: usize) {
    if size == 0 {
        panic!("Error: The container is empty.");
    }
}

/// Check that two dimensions are equal.
///
/// # Panics
/// Panics if `s1 != s2`.
#[inline]
pub fn check_size(s1: usize, s2: usize) {
    if s1 != s2 {
        panic!("Error: The dimensions mismatch.");
    }
}

/// Check if the matrix is a square matrix.
///
/// # Panics
/// Panics if the matrix is not square.
#[inline]
pub fn check_square(m: &Matrix) {
    if m.row_size() != m.col_size() {
        panic!("Error: The matrix is not a square matrix.");
    }
}
