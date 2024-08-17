use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use crate::{utility, Vector};

use pyinrs::Fraction;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Matrix {
    rows: Vec<Vector>,
}

impl Matrix {
    /// Create a new matrix object.
    pub fn new() -> Matrix {
        Matrix { rows: Vec::new() }
    }

    /// Create a row x col matrix with all identical elements.
    pub fn create(row: usize, col: usize, value: Fraction) -> Matrix {
        let mut rows = Vec::with_capacity(row);
        for _ in 0..row {
            rows.push(Vector::create(col, value));
        }
        Matrix { rows }
    }

    /// Create a row x col matrix with all 0 elements.
    pub fn zeros(row: usize, col: usize) -> Matrix {
        Matrix::create(row, col, 0.into())
    }

    /// Create a row x col matrix with all 1 elements.
    pub fn ones(row: usize, col: usize) -> Matrix {
        Matrix::create(row, col, 1.into())
    }

    /// Generate an n-order unit matrix.
    pub fn eye(n: usize) -> Matrix {
        let mut m = Matrix::zeros(n, n);
        for i in 0..n {
            m[i][i] = 1.into();
        }
        m
    }

    /// Return the number of rows in the matrix.
    pub fn row_size(&self) -> usize {
        self.rows.len()
    }

    /// Return the number of columns in the matrix.
    pub fn col_size(&self) -> usize {
        if self.row_size() == 0 {
            0
        } else {
            self.rows[0].size()
        }
    }

    /// Returns `true` if the matrix contains no elements.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Returns an iterator over the matrix.
    pub fn iter(&self) -> std::slice::Iter<Vector> {
        self.rows.iter()
    }

    /// Split this matrix by rows.
    pub fn split_row(&self, n: usize) -> (Matrix, Matrix) {
        utility::check_bounds(n, 0, self.row_size());

        let (mut first, mut second) = (Matrix::new(), Matrix::new());
        first.rows = self.rows[0..n].to_vec();
        second.rows = self.rows[n..].to_vec();

        (first, second)
    }

    /// Split this matrix by columns.
    pub fn split_col(&self, n: usize) -> (Matrix, Matrix) {
        utility::check_bounds(n, 0, self.col_size());

        let (mut first, mut second) = (Matrix::new(), Matrix::new());
        first.rows.resize(self.row_size(), Default::default());
        second.rows.resize(self.row_size(), Default::default());
        for r in 0..self.row_size() {
            first.rows[r].elements = self.rows[r].elements[..n].to_vec();
            second.rows[r].elements = self.rows[r].elements[n..].to_vec();
        }

        (first, second)
    }

    /// Returns the transpose of the matrix.
    pub fn transpose(&self) -> Matrix {
        let mut result = Matrix::zeros(self.col_size(), self.row_size());

        for i in 0..self.row_size() {
            for j in 0..self.col_size() {
                result[j][i] = self[i][j];
            }
        }
        result
    }

    /// Calculate the rank of this matrix.
    pub fn rank(&self) -> usize {
        let zeros = self.clone().to_row_echelon().rows.iter().filter(|row| row.is_zero()).count();
        self.row_size() - zeros
    }

    /// Calculate the determinant of this matrix.
    pub fn det(&self) -> Option<Fraction> {
        // check square matrix
        if self.row_size() != self.col_size() {
            return None;
        }

        let mut echelon = self.clone();
        echelon.to_row_echelon();
        let mut determinant = Fraction::from(1);
        for i in 0..echelon.row_size() {
            determinant *= echelon[i][i];
        }
        Some(determinant)
    }

    /// Calculate the inverse of this matrix.
    pub fn inv(&self) -> Option<Matrix> {
        // check square matrix and check invertible matrix
        if self.row_size() != self.col_size() || self.rank() != self.row_size() {
            return None;
        }

        let mut echelon = self.clone();
        // generate augmented matrix [A:E]
        echelon.expand_col(Matrix::eye(self.row_size()));
        // transforming [A:E] into a row echelon matrix
        echelon.to_row_echelon();
        // transform A into a diagonal matrix
        for c in 0..echelon.row_size() {
            for r in 0..c {
                echelon.e_row_sum(r, c, -(echelon[r][c] / echelon[c][c]));
            }
        }
        // transform A into a unit matrix
        for r in 0..echelon.row_size() {
            echelon.e_scalar_multiplication(r, Fraction::from(1) / echelon[r][r]);
        }
        // at this point, the original E is the inverse of A
        Some(echelon.split_col(self.row_size()).1)
    }

    /// Return the matrix that removed the i-th row and j-th column, 0 <= i, j < n.
    pub fn submatrix(&self, i: usize, j: usize) -> Matrix {
        let mut submatrix = Vec::with_capacity(self.row_size() - 1);
        for r in 0..self.row_size() {
            if r != i {
                let mut row = Vec::with_capacity(self.col_size() - 1);
                row.extend_from_slice(&self[r].elements[..j]);
                row.extend_from_slice(&self[r].elements[j + 1..]);
                submatrix.push(row);
            }
        }
        Matrix::from(submatrix)
    }

    /// Return the minor matrix.
    pub fn minor(&self) -> Matrix {
        let mut m = Matrix::zeros(self.row_size(), self.col_size());
        for r in 0..m.row_size() {
            for c in 0..m.col_size() {
                m[r][c] = self.submatrix(r, c).det().unwrap();
            }
        }
        m
    }

    /// Return the co-factor matrix.
    pub fn cofactor(&self) -> Matrix {
        let mut m = self.minor();
        for r in 0..m.row_size() {
            for c in 0..self.col_size() {
                // a11 -> a00, r+c parity unchanged
                if (r + c) & 1 == 1 {
                    m[r][c] = -m[r][c];
                }
            }
        }
        m
    }

    /// Return the adjugate matrix.
    pub fn adj(&self) -> Matrix {
        self.cofactor().transpose()
    }

    /// Returns `true` if the matrix is symmetric.
    pub fn is_symmetric(&self) -> bool {
        if self.row_size() != self.col_size() {
            return false;
        }

        for r in 0..self.row_size() {
            for c in 0..r {
                if self.rows[r][c] != self.rows[c][r] {
                    return false;
                }
            }
        }

        true
    }

    /// Expand this matrix by rows.
    pub fn expand_row(&mut self, mut matrix: Matrix) -> &Self {
        utility::check_size(self.col_size(), matrix.col_size());

        self.rows.append(&mut matrix.rows);
        self
    }

    /// Expand this matrix by columns.
    pub fn expand_col(&mut self, mut matrix: Matrix) -> &Self {
        utility::check_size(self.row_size(), matrix.row_size());

        for i in 0..self.row_size() {
            self.rows[i].elements.append(&mut matrix[i].elements);
        }
        self
    }

    /// Elementary Row Operations: Row Swap.
    pub fn e_row_swap(&mut self, i: usize, j: usize) -> &Self {
        self.rows.swap(i, j);
        self
    }

    /// Elementary Row Operations: Scalar Multiplication.
    pub fn e_scalar_multiplication(&mut self, i: usize, k: Fraction) -> &Self {
        self.rows[i] *= k;
        self
    }

    /// Elementary Row Operations: Row Sum.
    pub fn e_row_sum(&mut self, i: usize, j: usize, k: Fraction) -> &Self {
        let row = self[j].clone();
        self.rows[i] += row * k;
        self
    }

    /// Transform this matrix to general row echelon form.
    pub fn to_row_echelon(&mut self) -> &Self {
        // step 1: Gaussian elimination
        for i in 0..self.row_size() {
            let mut j: usize = 0;
            while j < self.col_size() && self.rows[i][j] == 0.into() {
                j += 1;
            }
            for k in i + 1..self.row_size() {
                if j < self.col_size() && self.rows[i][j] != 0.into() {
                    self.e_row_sum(k, i, -(self.rows[k][j] / self.rows[i][j]));
                }
            }
        }

        // step 2: transform to the row echelon form. It's so elegant, I'm a genius haha.
        self.rows.sort_by_key(|r| r.count_leading_zeros());
        self
    }
}

impl<const R: usize, const C: usize> From<[[Fraction; C]; R]> for Matrix {
    fn from(value: [[Fraction; C]; R]) -> Self {
        let rows = Vec::from(value.map(Vector::from));
        Matrix { rows }
    }
}

impl<const R: usize, const C: usize> From<[[i32; C]; R]> for Matrix {
    fn from(value: [[i32; C]; R]) -> Self {
        let rows = Vec::from(value.map(Vector::from));
        Matrix { rows }
    }
}

impl From<Vec<Vec<Fraction>>> for Matrix {
    fn from(value: Vec<Vec<Fraction>>) -> Self {
        let rows = value.into_iter().map(Vector::from).collect();
        Matrix { rows }
    }
}

impl From<Vec<Vec<i32>>> for Matrix {
    fn from(value: Vec<Vec<i32>>) -> Self {
        let rows = value.into_iter().map(Vector::from).collect();
        Matrix { rows }
    }
}

impl Index<usize> for Matrix {
    type Output = Vector;

    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}

impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.rows[index]
    }
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "[")?;

        // calc the max width of element
        let mut width = 0;
        for i in 0..self.row_size() {
            for j in 0..self.col_size() {
                width = width.max(format!("{}", self[i][j]).len());
            }
        }

        // align right, fill with space
        for i in 0..self.row_size() {
            for j in 0..self.col_size() {
                if j != 0 {
                    write!(f, " ")?;
                }
                write!(f, "{:>width$}", format!("{}", self[i][j]))?;
            }
            writeln!(f)?;
        }

        write!(f, "]")
    }
}

auto_ops::impl_op_ex!(+=|a: &mut Matrix, b: &Matrix| {
    utility::check_size(a.row_size(), b.row_size());
    utility::check_size(a.col_size(), b.col_size());

    for r in 0..a.row_size() {
        a[r] += &b[r];
    }
});

auto_ops::impl_op_ex!(+|a: &Matrix, b: &Matrix| -> Matrix {
    let mut a = a.clone();
    a += b;
    a
});

auto_ops::impl_op_ex!(-=|a: &mut Matrix, b: &Matrix| {
    utility::check_size(a.row_size(), b.row_size());
    utility::check_size(a.col_size(), b.col_size());

    for r in 0..a.row_size() {
        a[r] -= &b[r];
    }
});

auto_ops::impl_op_ex!(-|a: &Matrix, b: &Matrix| -> Matrix {
    let mut a = a.clone();
    a -= b;
    a
});

auto_ops::impl_op_ex!(*=|a: &mut Matrix, b: Fraction| {
    for r in 0..a.row_size() {
        a.rows[r] *= b;
    }
});

auto_ops::impl_op_ex_commutative!(*|a: Matrix, b: Fraction| -> Matrix {
    let mut a = a;
    a *= b;
    a
});

auto_ops::impl_op_ex!(*|a: &Matrix, b: &Matrix| -> Matrix {
    utility::check_size(a.col_size(), b.row_size());

    let mut result = Matrix::zeros(a.row_size(), b.col_size());
    let rt = b.transpose();
    for r in 0..a.row_size() {
        for c in 0..b.col_size() {
            result[r][c] = &a[r] * &rt[c];
        }
    }
    result
});
