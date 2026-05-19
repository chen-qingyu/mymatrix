use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use crate::{detail, MatrixError, Vector};

use pyinrs::Fraction;

/// Matrix with fractions as elements.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Matrix {
    rows: Vec<Vector>,
}

impl Matrix {
    /// Create a new empty matrix (0 rows, 0 columns).
    pub fn new() -> Self {
        Self { rows: Vec::new() }
    }

    /// Create a `row x col` matrix filled with `value`.
    pub fn create(row: usize, col: usize, value: Fraction) -> Self {
        let mut rows = Vec::with_capacity(row);
        for _ in 0..row {
            rows.push(Vector::create(col, value));
        }
        Self { rows }
    }

    /// Create a `row x col` zero matrix.
    pub fn zeros(row: usize, col: usize) -> Self {
        Self::create(row, col, 0.into())
    }

    /// Create a `row x col` matrix filled with 1.
    pub fn ones(row: usize, col: usize) -> Self {
        Self::create(row, col, 1.into())
    }

    /// Generate an `n x n` identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
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

    /// Returns `true` if the matrix is square.
    pub fn is_square(&self) -> bool {
        self.row_size() == self.col_size()
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

    /// Check if the matrix is upper triangular matrix.
    pub fn is_upper(&self) -> bool {
        if self.row_size() != self.col_size() {
            return false;
        }

        for r in 1..self.row_size() {
            for c in 0..r {
                if self[r][c] != 0.into() {
                    return false;
                }
            }
        }

        true
    }

    /// Check if the matrix is lower triangular matrix.
    pub fn is_lower(&self) -> bool {
        if self.row_size() != self.col_size() {
            return false;
        }

        for c in 1..self.col_size() {
            for r in 0..c {
                if self[r][c] != 0.into() {
                    return false;
                }
            }
        }

        true
    }

    /// Check if the matrix is diagonal matrix.
    pub fn is_diagonal(&self) -> bool {
        self.is_lower() && self.is_upper()
    }

    /// Calculate the trace of the matrix.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn trace(&self) -> Fraction {
        detail::check_square(self);

        let mut tr = Fraction::new();
        for i in 0..self.row_size() {
            tr += self[i][i];
        }

        tr
    }

    /// Returns the transpose of the matrix.
    pub fn transpose(&self) -> Self {
        let mut result = Self::zeros(self.col_size(), self.row_size());

        for i in 0..self.row_size() {
            for j in 0..self.col_size() {
                result[j][i] = self[i][j];
            }
        }
        result
    }

    /// Transform this matrix to general row echelon form.
    pub fn row_echelon_form(&self) -> Self {
        let mut m = self.clone();

        // Gaussian elimination
        for i in 0..m.row_size() {
            let mut j: usize = 0;
            while j < m.col_size() && m.rows[i][j] == 0.into() {
                j += 1;
            }
            for k in i + 1..m.row_size() {
                if j < m.col_size() && m.rows[i][j] != 0.into() {
                    m.e_row_sum(k, i, -m.rows[k][j] / m.rows[i][j]);
                }
            }
        }

        // transform to the row echelon form. It's so elegant, I'm a genius haha.
        m.rows.sort_by_key(|r| r.count_leading_zeros());

        m
    }

    /// Transform this matrix to reduced row echelon form.
    pub fn row_canonical_form(&self) -> Self {
        let mut m = self.row_echelon_form();

        // find the pivot column for each row
        let pivot_cols: Vec<Option<usize>> = (0..m.row_size())
            .map(|r| if m[r].is_zero() { None } else { Some(m[r].count_leading_zeros()) })
            .collect();

        // eliminate elements above each pivot
        for (pivot_row, &pivot_col) in pivot_cols.iter().enumerate() {
            if let Some(col) = pivot_col {
                for r in 0..pivot_row {
                    if m[r][col] != 0.into() {
                        m.e_row_sum(r, pivot_row, -(m[r][col] / m[pivot_row][col]));
                    }
                }
            }
        }

        // make each pivot equal to 1
        for r in 0..m.row_size() {
            if let Some(col) = pivot_cols[r] {
                if m[r][col] != 0.into() {
                    m.e_scalar_multiplication(r, Fraction::from(1) / m[r][col]);
                }
            }
        }

        m
    }

    /// Calculate the determinant of this matrix.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn det(&self) -> Fraction {
        detail::check_square(self);

        let n = self.row_size();
        let mut a = self.clone();

        let mut det = Fraction::from(1);
        for i in 0..n {
            let mut pivot = i;
            for j in i + 1..n {
                if a[j][i].abs() > a[pivot][i].abs() {
                    pivot = j;
                }
            }
            if pivot != i {
                a.e_row_swap(i, pivot);
                det = -det;
            }
            if a[i][i] == 0.into() {
                return Fraction::new();
            }
            det *= a[i][i];
            for j in i + 1..n {
                a.e_row_sum(j, i, -a[j][i] / a[i][i]);
            }
        }
        det
    }

    /// Return the matrix obtained by removing the `i`-th row and `j`-th column.
    ///
    /// # Panics
    /// Panics if `i` or `j` is out of bounds.
    pub fn submatrix(&self, i: usize, j: usize) -> Self {
        let mut submatrix = Vec::with_capacity(self.row_size() - 1);
        for r in 0..self.row_size() {
            if r != i {
                let mut row = Vec::with_capacity(self.col_size() - 1);
                row.extend_from_slice(&self[r].elements[..j]);
                row.extend_from_slice(&self[r].elements[j + 1..]);
                submatrix.push(row);
            }
        }
        Self::from(submatrix)
    }

    /// Return the minor matrix.
    pub fn minor(&self) -> Self {
        let mut m = Self::zeros(self.row_size(), self.col_size());
        for r in 0..m.row_size() {
            for c in 0..m.col_size() {
                m[r][c] = self.submatrix(r, c).det();
            }
        }
        m
    }

    /// Return the cofactor matrix.
    pub fn cofactor(&self) -> Self {
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
    pub fn adj(&self) -> Self {
        self.cofactor().transpose()
    }

    /// Calculate the inverse of this matrix.
    ///
    /// Returns `Err(MatrixError::Singular)` if the matrix is not invertible.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn inv(&self) -> Result<Self, MatrixError> {
        detail::check_square(self);

        // inverse of empty matrix is empty matrix
        if self.is_empty() {
            return Ok(Matrix::new());
        }

        // generate augmented matrix [A:E] and transform [A:E] to reduced row echelon form and split
        let n = self.row_size();
        let rref = self.clone().expand_col(Self::identity(n)).row_canonical_form().split_col(n);

        // now, the original E is the inverse of A if rank = n
        if !rref.0[n - 1].is_zero() {
            Ok(rref.1)
        } else {
            Err(MatrixError::Singular)
        }
    }

    /// Calculate the rank of this matrix.
    pub fn rank(&self) -> usize {
        let zeros = self.row_echelon_form().rows.iter().filter(|row| row.is_zero()).count();
        self.row_size() - zeros
    }

    /// LU decomposition using the Doolittle algorithm.
    ///
    /// Returns `Err(MatrixError::Singular)` if the matrix is singular.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn lu_decomposition(&self) -> Result<(Self, Self), MatrixError> {
        detail::check_square(self);

        let n = self.row_size();

        if self.is_upper() {
            // include zero matrix
            return Ok((Matrix::identity(n), self.clone()));
        }

        let mut l = Self::identity(n);
        let mut u = Self::zeros(n, n);

        for i in 0..n {
            for j in 0..(i + 1) {
                let mut sum = Fraction::new();
                for k in 0..j {
                    sum += l[j][k] * u[k][i];
                }
                u[j][i] = self[j][i] - sum;
            }

            for j in (i + 1)..n {
                let mut sum = Fraction::new();
                for k in 0..i {
                    sum += l[j][k] * u[k][i];
                }
                if u[i][i] == 0.into() {
                    return Err(MatrixError::Singular);
                }
                l[j][i] = (self[j][i] - sum) / u[i][i];
            }
        }

        Ok((l, u))
    }

    /// Calculate the integer power of a square matrix using binary exponentiation.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn pow(&self, exp: u32) -> Self {
        detail::check_square(self);

        match exp {
            0 => Self::identity(self.row_size()),
            1 => self.clone(),
            _ => {
                let mut result = Self::identity(self.row_size());
                let mut base = self.clone();
                let mut e = exp;
                while e > 0 {
                    if e & 1 == 1 {
                        result = &result * &base;
                    }
                    base = &base * &base;
                    e >>= 1;
                }
                result
            }
        }
    }

    /// Solve the linear system `Ax = b`, where `A` is this square matrix.
    ///
    /// Returns `Err(MatrixError::Singular)` if `A` is singular
    /// (no unique solution or inconsistent system).
    ///
    /// # Panics
    /// Panics if the matrix is not square, or if the size of `b` does not match
    /// the number of rows of `A`.
    pub fn solve(&self, b: &Vector) -> Result<Vector, MatrixError> {
        detail::check_square(self);
        detail::check_size(self.row_size(), b.size());

        // build augmented matrix [A | b] and compute RREF
        let mut augmented = self.clone();
        for r in 0..self.row_size() {
            augmented.rows[r].elements.push(b[r]);
        }
        let rref = augmented.row_canonical_form();

        // check for inconsistency: row of [0 ... 0 | non-zero]
        for r in 0..rref.row_size() {
            let all_zero = rref[r].elements[..rref.col_size() - 1].iter().all(|&x| x == 0.into());
            if all_zero && rref[r][rref.col_size() - 1] != 0.into() {
                return Err(MatrixError::Singular);
            }
        }

        // extract solution (works when rank = n)
        let n = self.row_size();
        let mut x = Vector::zeros(n);
        for i in 0..n {
            // check if column i is a pivot column
            if rref[i][i] != 0.into() {
                x[i] = rref[i][n] / rref[i][i];
            } else {
                // free variable; no unique solution
                return Err(MatrixError::Singular);
            }
        }
        Ok(x)
    }

    /// Split this matrix by rows.
    pub fn split_row(&self, n: usize) -> (Self, Self) {
        detail::check_bounds(n, 0, self.row_size());

        let (mut first, mut second) = (Self::new(), Self::new());
        first.rows = self.rows[0..n].to_vec();
        second.rows = self.rows[n..].to_vec();

        (first, second)
    }

    /// Split this matrix by columns.
    pub fn split_col(&self, n: usize) -> (Self, Self) {
        detail::check_bounds(n, 0, self.col_size());

        let (mut first, mut second) = (Self::new(), Self::new());
        first.rows.resize(self.row_size(), Default::default());
        second.rows.resize(self.row_size(), Default::default());
        for r in 0..self.row_size() {
            first.rows[r].elements = self.rows[r].elements[..n].to_vec();
            second.rows[r].elements = self.rows[r].elements[n..].to_vec();
        }

        (first, second)
    }

    /// Expand this matrix by rows.
    pub fn expand_row(&mut self, mut matrix: Self) -> &Self {
        detail::check_size(self.col_size(), matrix.col_size());

        self.rows.append(&mut matrix.rows);
        self
    }

    /// Expand this matrix by columns.
    pub fn expand_col(&mut self, mut matrix: Self) -> &Self {
        detail::check_size(self.row_size(), matrix.row_size());

        for i in 0..self.row_size() {
            self.rows[i].elements.append(&mut matrix[i].elements);
        }
        self
    }

    /// Elementary Row Operations: Row Swap. (A[i] <=> A[j])
    pub fn e_row_swap(&mut self, i: usize, j: usize) -> &Self {
        self.rows.swap(i, j);
        self
    }

    /// Elementary Row Operations: Scalar Multiplication. (A[i] *= k)
    pub fn e_scalar_multiplication(&mut self, i: usize, k: Fraction) -> &Self {
        self.rows[i] *= k;
        self
    }

    /// Elementary Row Operations: Row Sum. (A[i] += A[j] * k)
    pub fn e_row_sum(&mut self, i: usize, j: usize, k: Fraction) -> &Self {
        // clone row j to avoid aliasing with self modification, then inline scalar multiply
        // to avoid a second clone inside the `*` operator
        let scaled = {
            let mut row = self[j].clone();
            for elem in &mut row.elements {
                *elem *= k;
            }
            row
        };
        self.rows[i] += scaled;
        self
    }
}

impl<const R: usize, const C: usize> From<[[Fraction; C]; R]> for Matrix {
    fn from(value: [[Fraction; C]; R]) -> Self {
        let rows = Vec::from(value.map(Vector::from));
        Self { rows }
    }
}

impl<const R: usize, const C: usize> From<[[i32; C]; R]> for Matrix {
    fn from(value: [[i32; C]; R]) -> Self {
        let rows = Vec::from(value.map(Vector::from));
        Self { rows }
    }
}

impl From<Vec<Vec<Fraction>>> for Matrix {
    fn from(value: Vec<Vec<Fraction>>) -> Self {
        let rows = value.into_iter().map(Vector::from).collect();
        Self { rows }
    }
}

impl From<Vec<Vec<i32>>> for Matrix {
    fn from(value: Vec<Vec<i32>>) -> Self {
        let rows = value.into_iter().map(Vector::from).collect();
        Self { rows }
    }
}

impl From<Vec<Vector>> for Matrix {
    fn from(value: Vec<Vector>) -> Self {
        Self { rows: value }
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
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
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
    detail::check_size(a.row_size(), b.row_size());
    detail::check_size(a.col_size(), b.col_size());

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
    detail::check_size(a.row_size(), b.row_size());
    detail::check_size(a.col_size(), b.col_size());

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

auto_ops::impl_op_ex!(*=|a: &mut Matrix, b: i32| {
    for r in 0..a.row_size() {
        a.rows[r] *= b;
    }
});

auto_ops::impl_op_ex_commutative!(*|a: Matrix, b: i32| -> Matrix {
    let mut a = a;
    a *= b;
    a
});

auto_ops::impl_op_ex!(/=|a: &mut Matrix, b: Fraction| {
    for r in 0..a.row_size() {
        a.rows[r] /= b;
    }
});

auto_ops::impl_op_ex!(/|a: &Matrix, b: Fraction| -> Matrix {
    let mut a = a.clone();
    a /= b;
    a
});

auto_ops::impl_op_ex!(*|a: &Matrix, b: &Matrix| -> Matrix {
    detail::check_size(a.col_size(), b.row_size());

    let mut result = Matrix::zeros(a.row_size(), b.col_size());
    let rt = b.transpose();
    for r in 0..a.row_size() {
        for c in 0..b.col_size() {
            result[r][c] = &a[r] * &rt[c];
        }
    }
    result
});

auto_ops::impl_op_ex!(*|a: &Matrix, b: &Vector| -> Vector {
    detail::check_size(a.col_size(), b.size());

    let mut result = Vector::zeros(a.row_size());
    for r in 0..a.row_size() {
        result[r] = &a[r] * b;
    }
    result
});

impl std::ops::Neg for Matrix {
    type Output = Matrix;

    fn neg(mut self) -> Matrix {
        for row in &mut self.rows {
            for elem in &mut row.elements {
                *elem = -(*elem);
            }
        }
        self
    }
}

impl std::ops::Neg for &Matrix {
    type Output = Matrix;

    fn neg(self) -> Matrix {
        -(self.clone())
    }
}

impl IntoIterator for Matrix {
    type Item = Vector;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.rows.into_iter()
    }
}
