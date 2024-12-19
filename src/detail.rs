use crate::Matrix;

// Check whether the index is valid (begin <= pos < end).
#[inline]
pub fn check_bounds(pos: usize, begin: usize, end: usize) {
    if pos < begin || pos >= end {
        panic!("Error: Index out of range.");
    }
}

// Check whether is not empty.
#[inline]
pub fn check_empty(size: usize) {
    if size == 0 {
        panic!("Error: The container is empty.");
    }
}

// Check that two vectors are of the same size.
#[inline]
pub fn check_size(s1: usize, s2: usize) {
    if s1 != s2 {
        panic!("Error: The dimensions mismatch.");
    }
}

// Check if the matrix is a square matrix.
#[inline]
pub fn check_square(m: &Matrix) {
    if m.row_size() != m.col_size() {
        panic!("Error: This function applies only to square matrices.");
    }
}
