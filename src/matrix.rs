use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use crate::{detail, Vector};

use pyinrs::Fraction;

/// Matrix with fractions as elements.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Matrix {
    rows: Vec<Vector>,
}

impl Matrix {
    /// Create a new matrix object.
    pub fn new() -> Self {
        Self { rows: Vec::new() }
    }

    /// Create a row x col matrix with all identical elements.
    pub fn create(row: usize, col: usize, value: Fraction) -> Self {
        let mut rows = Vec::with_capacity(row);
        for _ in 0..row {
            rows.push(Vector::create(col, value));
        }
        Self { rows }
    }

    /// Create a row x col matrix with all 0 elements.
    pub fn zeros(row: usize, col: usize) -> Self {
        Self::create(row, col, 0.into())
    }

    /// Create a row x col matrix with all 1 elements.
    pub fn ones(row: usize, col: usize) -> Self {
        Self::create(row, col, 1.into())
    }

    /// Generate an n-order identity matrix.
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

        return true;
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

        return true;
    }

    /// Check if the matrix is diagonal matrix.
    pub fn is_diagonal(&self) -> bool {
        return self.is_lower() && self.is_upper();
    }

    /// Calculate the trace of the matrix.
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

        // step 1: Gaussian elimination
        for i in 0..m.row_size() {
            let mut j: usize = 0;
            while j < m.col_size() && m.rows[i][j] == 0.into() {
                j += 1;
            }
            for k in i + 1..m.row_size() {
                if j < m.col_size() && m.rows[i][j] != 0.into() {
                    m.e_row_sum(k, i, -(m.rows[k][j] / m.rows[i][j]));
                }
            }
        }

        // step 2: transform to the row echelon form. It's so elegant, I'm a genius haha.
        m.rows.sort_by_key(|r| r.count_leading_zeros());

        m
    }

    /// Transform this matrix to reduced row echelon form.
    pub fn row_canonical_form(&self) -> Self {
        let mut m = self.row_echelon_form();

        let n = usize::min(m.row_size(), m.col_size());

        // step 3: eliminate elements above the pivot
        for c in 0..n {
            for r in 0..c {
                if m[c][c] != 0.into() {
                    m.e_row_sum(r, c, -(m[r][c] / m[c][c]));
                }
            }
        }

        // step 4: make pivot equals 1
        for r in 0..n {
            if m[r][r] != 0.into() {
                m.e_scalar_multiplication(r, Fraction::from(1) / m[r][r]);
            }
        }

        m
    }

    /// Calculate the determinant of this matrix.
    pub fn det(&self) -> Fraction {
        detail::check_square(self);

        let u = self.row_echelon_form();
        let mut determinant = Fraction::from(1);
        for i in 0..u.row_size() {
            determinant *= u[i][i];
        }
        determinant
    }

    /// Return the matrix that removed the i-th row and j-th column, 0 <= i, j < n.
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
    pub fn inv(&self) -> Option<Self> {
        detail::check_square(self);

        // check empty
        if self.is_empty() {
            return Some(Matrix::new());
        }

        // check invertible
        if self.rank() != self.row_size() {
            return None;
        }

        // generate augmented matrix [A:E] and transform [A:E] to reduced row echelon form
        let echelon = self.clone().expand_col(Self::identity(self.row_size())).row_canonical_form();

        // at this point, the original E is the inverse of A
        Some(echelon.split_col(self.row_size()).1)
    }

    /// Calculate the rank of this matrix.
    pub fn rank(&self) -> usize {
        let zeros = self.row_echelon_form().rows.iter().filter(|row| row.is_zero()).count();
        self.row_size() - zeros
    }

    /// LU decomposition, use Doolittle algorithm.
    pub fn lu_decomposition(&self) -> (Self, Self) {
        detail::check_square(self);

        let n = self.row_size();

        if self.is_upper() {
            return (Matrix::zeros(n, n), self.clone());
        } else if self.is_lower() {
            return (self.clone(), Matrix::zeros(n, n));
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
                l[j][i] = (self[j][i] - sum) / u[i][i];
            }
        }

        (l, u)
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

impl IntoIterator for Matrix {
    type Item = Vector;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.rows.into_iter()
    }
}
