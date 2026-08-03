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
        Self {
            rows: vec![Vector::create(col, value); row],
        }
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

    /// Extract the main diagonal as a vector.
    ///
    /// For a non-square matrix, returns `min(row_size, col_size)` elements.
    pub fn diag(&self) -> Vector {
        let n = self.row_size().min(self.col_size());
        Vector::from((0..n).map(|i| self[i][i]).collect::<Vec<_>>())
    }

    /// Create a diagonal matrix from the elements of a vector.
    pub fn from_diagonal(v: &Vector) -> Self {
        let mut m = Self::zeros(v.size(), v.size());
        for i in 0..v.size() {
            m[i][i] = v[i];
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

    /// Return the `c`-th column as a vector.
    ///
    /// # Panics
    /// Panics if `c` is out of bounds.
    pub fn col(&self, c: usize) -> Vector {
        Vector::from((0..self.row_size()).map(|r| self[r][c]).collect::<Vec<_>>())
    }

    /// Returns `true` if the matrix has no rows (i.e. is `0 x 0`).
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

    /// Rotate the matrix 90 degrees counter-clockwise.
    pub fn rotate_left(&self) -> Self {
        let m = self.col_size();
        let mut result = Self::zeros(m, self.row_size());

        for r in 0..self.row_size() {
            for c in 0..m {
                result[m - c - 1][r] = self[r][c];
            }
        }
        result
    }

    /// Rotate the matrix 90 degrees clockwise.
    pub fn rotate_right(&self) -> Self {
        let n = self.row_size();
        let mut result = Self::zeros(self.col_size(), n);

        for r in 0..n {
            for c in 0..self.col_size() {
                result[c][n - r - 1] = self[r][c];
            }
        }
        result
    }

    /// Calculate the Kronecker product of this matrix with `that`.
    pub fn kron(&self, that: &Self) -> Self {
        let (m, n) = (self.row_size(), self.col_size());
        let (p, q) = (that.row_size(), that.col_size());
        let mut result = Self::zeros(m * p, n * q);
        for i in 0..m {
            for j in 0..n {
                for k in 0..p {
                    for l in 0..q {
                        result[i * p + k][j * q + l] = self[i][j] * that[k][l];
                    }
                }
            }
        }
        result
    }

    /// Calculate the Hadamard (element-wise) product of this matrix with `that`.
    ///
    /// # Panics
    /// Panics if the dimensions do not match.
    pub fn hadamard(&self, that: &Self) -> Self {
        detail::check_size(self.row_size(), that.row_size());
        detail::check_size(self.col_size(), that.col_size());

        let mut result = self.clone();
        for r in 0..result.row_size() {
            for c in 0..result.col_size() {
                result[r][c] *= that[r][c];
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
            if j < m.col_size() {
                for k in i + 1..m.row_size() {
                    m.e_row_sum(k, i, -m.rows[k][j] / m.rows[i][j]);
                }
            }
        }

        // order rows by pivot column; zero rows (key = size) sink to the bottom
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
                m.e_scalar_multiplication(r, Fraction::from(1) / m[r][col]);
            }
        }

        m
    }

    /// Calculate the determinant of this matrix.
    ///
    /// The determinant of the empty (`0 x 0`) matrix is `1` by convention.
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
    /// Panics if `i` or `j` is out of bounds, or if the matrix has no row or
    /// no column to remove (including an empty matrix).
    pub fn submatrix(&self, i: usize, j: usize) -> Self {
        // an empty (or row-less/column-less) matrix has nothing to remove,
        // so any index is out of range; avoids underflowing `size - 1`
        if self.row_size() == 0 || self.col_size() == 0 {
            panic!("Error: Index out of range.");
        }
        detail::check_bounds(i, 0, self.row_size() - 1);
        detail::check_bounds(j, 0, self.col_size() - 1);

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
    ///
    /// # Panics
    /// Panics if the matrix is not square.
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
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn cofactor(&self) -> Self {
        let mut m = self.minor();
        for r in 0..m.row_size() {
            for c in 0..self.col_size() {
                // sign of each cofactor is (-1)^(r+c)
                if (r + c) & 1 == 1 {
                    m[r][c] = -m[r][c];
                }
            }
        }
        m
    }

    /// Return the adjugate matrix.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn adj(&self) -> Self {
        self.cofactor().transpose()
    }

    /// Calculate the inverse of this matrix.
    ///
    /// # Errors
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

    /// Calculate the Moore–Penrose pseudo-inverse of this matrix.
    ///
    /// Generalizes `inv` to non-square or singular matrices. Uses the exact
    /// rank factorization `A = B * C` and the formula
    /// `A⁺ = Cᵀ (C Cᵀ)⁻¹ (Bᵀ B)⁻¹ Bᵀ`.
    pub fn pseudo_inverse(&self) -> Self {
        let b = self.col_space().transpose(); // m x r
        if b.col_size() == 0 {
            return Self::zeros(self.col_size(), self.row_size());
        }

        let bt_b_inv = (&b.transpose() * &b).inv().unwrap();
        let c = &bt_b_inv * &b.transpose() * self; // r x n
        let c_ct_inv = (&c * &c.transpose()).inv().unwrap();
        &c.transpose() * &c_ct_inv * &bt_b_inv * &b.transpose()
    }

    /// Calculate the rank of this matrix.
    pub fn rank(&self) -> usize {
        let zeros = self.row_echelon_form().rows.iter().filter(|row| row.is_zero()).count();
        self.row_size() - zeros
    }

    /// Return a basis of the row space, as rows of a matrix.
    pub fn row_space(&self) -> Self {
        let rows = self.row_echelon_form().rows.into_iter().filter(|r| !r.is_zero()).collect();
        Self { rows }
    }

    /// Return a basis of the column space, as rows of a matrix.
    ///
    /// The basis consists of the columns of this matrix at the pivot
    /// columns of its reduced row echelon form.
    pub fn col_space(&self) -> Self {
        let (_, pivot_cols) = self.rref_and_pivot_cols();
        let mut basis = Self::zeros(pivot_cols.len(), self.row_size());
        for (b, &pivot) in pivot_cols.iter().enumerate() {
            for r in 0..self.row_size() {
                basis[b][r] = self[r][pivot];
            }
        }
        basis
    }

    /// Return a basis of the null space, as rows of a matrix.
    ///
    /// Each row is a null vector; the result has `col_size()` columns and
    /// `col_size() - rank` rows.
    pub fn null_space(&self) -> Self {
        let (rref, pivot_cols) = self.rref_and_pivot_cols();
        Self::null_basis(&rref, &pivot_cols, self.col_size())
    }

    /// Rows of the null-space basis (one per free column) built from an RREF.
    fn null_basis(rref: &Self, pivot_cols: &[usize], col_count: usize) -> Self {
        let free_cols: Vec<usize> = (0..col_count).filter(|c| !pivot_cols.contains(c)).collect();
        let mut basis = Vec::with_capacity(free_cols.len());
        for &f in &free_cols {
            let mut v = Vector::zeros(col_count);
            v[f] = 1.into();
            for (i, &p) in pivot_cols.iter().enumerate() {
                v[p] = -rref[i][f];
            }
            basis.push(v);
        }
        Self::from(basis)
    }

    /// Return the reduced row echelon form and its pivot columns.
    fn rref_and_pivot_cols(&self) -> (Self, Vec<usize>) {
        let rref = self.row_canonical_form();
        let pivot_cols: Vec<usize> = (0..rref.row_size())
            .filter_map(|r| if rref[r].is_zero() { None } else { Some(rref[r].count_leading_zeros()) })
            .collect();
        (rref, pivot_cols)
    }

    /// LDL^T decomposition (rational Cholesky) for symmetric positive definite matrices.
    ///
    /// Returns `(L, d)`, where `L` is a unit lower triangular matrix and `d` is the
    /// diagonal of `D`, satisfying `A = L * D * L^T`.
    ///
    /// This is the exact-rational analogue of the classical Cholesky decomposition.
    /// Instead of $A = G G^T$ (which requires `sqrt`), it computes
    /// $A = L D L^T$ where $G = L \sqrt{D}$.
    ///
    /// # Errors
    /// - `Err(MatrixError::NotPositiveDefinite)` if the matrix is not symmetric
    ///   positive definite.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn cholesky(&self) -> Result<(Self, Vector), MatrixError> {
        detail::check_square(self);

        if !self.is_symmetric() {
            return Err(MatrixError::NotPositiveDefinite);
        }

        let n = self.row_size();

        // L is unit lower triangular (initialized as identity)
        let mut l = Self::identity(n);
        // Diagonal of D
        let mut d = Vector::zeros(n);

        for i in 0..n {
            // d[i] = A[i][i] - sum(L[i][k]^2 * d[k] for k in 0..i)
            let mut di = self[i][i];
            for k in 0..i {
                di -= l[i][k] * l[i][k] * d[k];
            }

            if di <= 0.into() {
                return Err(MatrixError::NotPositiveDefinite);
            }

            d[i] = di;

            // L[j][i] = (A[j][i] - sum(L[j][k]*L[i][k]*d[k] for k in 0..i)) / d[i]
            for j in (i + 1)..n {
                let mut lji = self[j][i];
                for k in 0..i {
                    lji -= l[j][k] * l[i][k] * d[k];
                }
                l[j][i] = lji / d[i];
            }
        }

        Ok((l, d))
    }

    /// LU decomposition using the Doolittle algorithm with row pivoting.
    ///
    /// Rows are swapped only when a pivot is zero, so `L * U` equals `A`
    /// after applying the same row swaps to `A`.
    ///
    /// # Errors
    /// Returns `Err(MatrixError::Singular)` if and only if the matrix is
    /// singular.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn lu_decomposition(&self) -> Result<(Self, Self), MatrixError> {
        detail::check_square(self);

        let n = self.row_size();
        let mut a = self.clone();

        for i in 0..n {
            // find a non-zero pivot in column i (rows i..n)
            let mut pivot = i;
            while pivot < n && a[pivot][i] == 0.into() {
                pivot += 1;
            }
            if pivot == n {
                return Err(MatrixError::Singular);
            }
            if pivot != i {
                a.rows.swap(i, pivot);
            }

            // column i of L (multipliers), then update the trailing submatrix (U)
            let pivot = a[i][i];
            for j in (i + 1)..n {
                a[j][i] /= pivot;
                for k in (i + 1)..n {
                    a[j][k] = a[j][k] - a[j][i] * a[i][k];
                }
            }
        }

        // extract L (unit lower triangular) and U (upper triangular)
        let mut l = Self::identity(n);
        let mut u = Self::zeros(n, n);
        for i in 0..n {
            for j in 0..i {
                l[i][j] = a[i][j];
            }
            for j in i..n {
                u[i][j] = a[i][j];
            }
        }

        Ok((l, u))
    }

    /// Calculate the non-negative integer power of a square matrix using
    /// binary exponentiation (`exp = 0` yields the identity matrix).
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn pow(&self, exp: u32) -> Self {
        detail::check_square(self);

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

    /// Calculate the characteristic polynomial of this matrix.
    ///
    /// Returns the coefficients `[c_0, c_1, ..., c_n]` of
    /// `det(λI - A) = c_0 λ^n + c_1 λ^(n-1) + ... + c_n`, with `c_0 = 1`.
    ///
    /// # Panics
    /// Panics if the matrix is not square.
    pub fn characteristic_polynomial(&self) -> Vec<Fraction> {
        detail::check_square(self);
        let n = self.row_size();

        // elementary symmetric polynomials of the eigenvalues, via Newton's identities
        let mut s = vec![Fraction::new(); n + 1]; // s[k] = trace(A^k)
        let mut e = vec![Fraction::new(); n + 1]; // e[k] = k-th elementary symmetric polynomial
        e[0] = 1.into();
        for k in 1..=n {
            s[k] = self.pow(k as u32).trace();
            let mut sum = Fraction::new();
            for i in 1..=k {
                let term = e[k - i] * s[i];
                sum += if i % 2 == 1 { term } else { -term };
            }
            e[k] = sum / Fraction::from(k as i32);
        }

        // coefficient of λ^(n-k) is (-1)^k e_k
        (0..=n).map(|k| if k % 2 == 0 { e[k] } else { -e[k] }).collect()
    }

    /// Return the general solution of the linear system `Ax = b`, where `A`
    /// is this square matrix.
    ///
    /// Returns `(x_p, N)`, where `x_p` is a particular solution and the rows
    /// of `N` form a basis of the null space, so every solution is
    /// `x_p + N[0]·c_1 + N[1]·c_2 + ...`. `N` is empty exactly when the
    /// solution is unique.
    ///
    /// # Errors
    /// Returns `Err(MatrixError::Inconsistent)` if the system is inconsistent
    /// (has no solution).
    ///
    /// # Panics
    /// Panics if the matrix is not square, or if the size of `b` does not match
    /// the number of rows of `A`.
    pub fn general_solution(&self, b: &Vector) -> Result<(Vector, Self), MatrixError> {
        detail::check_square(self);
        detail::check_size(self.row_size(), b.size());
        let n = self.row_size();

        // build augmented matrix [A | b] and compute RREF
        let mut augmented = self.clone();
        for r in 0..n {
            augmented.rows[r].elements.push(b[r]);
        }
        let rref = augmented.row_canonical_form();

        // check for inconsistency: row of [0 ... 0 | non-zero]
        for r in 0..rref.row_size() {
            let all_zero = rref[r].elements[..rref.col_size() - 1].iter().all(|&x| x == 0.into());
            if all_zero && rref[r][rref.col_size() - 1] != 0.into() {
                return Err(MatrixError::Inconsistent);
            }
        }

        // particular solution with free variables set to 0
        let pivot_cols: Vec<usize> = (0..rref.row_size())
            .filter_map(|r| if rref[r].is_zero() { None } else { Some(rref[r].count_leading_zeros()) })
            .collect();
        let mut xp = Vector::zeros(n);
        for (i, &p) in pivot_cols.iter().enumerate() {
            xp[p] = rref[i][n];
        }

        Ok((xp, Self::null_basis(&rref, &pivot_cols, n)))
    }

    /// Solve the linear system `Ax = b`, where `A` is this square matrix.
    ///
    /// Solve the linear system `Ax = b`, where `A` is this square matrix.
    ///
    /// Returns the unique solution.
    ///
    /// # Errors
    /// - `Err(MatrixError::Inconsistent)` if the system has no solution.
    /// - `Err(MatrixError::Singular)` if the system is consistent but has
    ///   infinitely many solutions (`A` is singular).
    ///
    /// For the general solution, use `general_solution`.
    ///
    /// # Panics
    /// Panics if the matrix is not square, or if the size of `b` does not match
    /// the number of rows of `A`.
    pub fn solve(&self, b: &Vector) -> Result<Vector, MatrixError> {
        let (xp, null) = self.general_solution(b)?;
        if null.is_empty() {
            Ok(xp)
        } else {
            Err(MatrixError::Singular)
        }
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

    /// Elementary Row Operations: Row Swap. (`A[i] <=> A[j]`)
    pub fn e_row_swap(&mut self, i: usize, j: usize) -> &Self {
        self.rows.swap(i, j);
        self
    }

    /// Elementary Row Operations: Scalar Multiplication. (`A[i] *= k`)
    pub fn e_scalar_multiplication(&mut self, i: usize, k: Fraction) -> &Self {
        self.rows[i] *= k;
        self
    }

    /// Elementary Row Operations: Row Sum. (`A[i] += A[j] * k`)
    pub fn e_row_sum(&mut self, i: usize, j: usize, k: Fraction) -> &Self {
        let scaled = self[j].clone() * k;
        self.rows[i] += scaled;
        self
    }

    /// Return an iterator over the rows.
    pub fn iter(&self) -> std::slice::Iter<'_, Vector> {
        self.rows.iter()
    }

    /// Return a mutable iterator over the rows.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Vector> {
        self.rows.iter_mut()
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
        Self::from(value.into_iter().map(Vector::from).collect::<Vec<_>>())
    }
}

impl From<Vec<Vec<i32>>> for Matrix {
    fn from(value: Vec<Vec<i32>>) -> Self {
        Self::from(value.into_iter().map(Vector::from).collect::<Vec<_>>())
    }
}

impl From<Vec<Vector>> for Matrix {
    fn from(value: Vec<Vector>) -> Self {
        if let Some(first) = value.first() {
            let len = first.size();
            for v in &value[1..] {
                assert_eq!(v.size(), len, "Error: All rows must have the same length.");
            }
        }
        Self { rows: value }
    }
}

impl FromIterator<Vector> for Matrix {
    fn from_iter<T: IntoIterator<Item = Vector>>(iter: T) -> Self {
        Self::from(iter.into_iter().collect::<Vec<_>>())
    }
}

impl FromIterator<Vec<Fraction>> for Matrix {
    fn from_iter<T: IntoIterator<Item = Vec<Fraction>>>(iter: T) -> Self {
        Self::from(iter.into_iter().collect::<Vec<_>>())
    }
}

impl FromIterator<Vec<i32>> for Matrix {
    fn from_iter<T: IntoIterator<Item = Vec<i32>>>(iter: T) -> Self {
        Self::from(iter.into_iter().collect::<Vec<_>>())
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
        use std::fmt::Write;

        writeln!(f, "[")?;

        // calc the max width of element
        let mut buf = String::new();
        let mut width = 0;
        for i in 0..self.row_size() {
            for j in 0..self.col_size() {
                write!(buf, "{}", self[i][j]).unwrap();
                width = width.max(buf.len());
                buf.clear();
            }
        }

        // align right, fill with space
        for i in 0..self.row_size() {
            for j in 0..self.col_size() {
                if j != 0 {
                    write!(f, " ")?;
                }
                write!(buf, "{}", self[i][j]).unwrap();
                write!(f, "{:>width$}", buf)?;
                buf.clear();
            }
            writeln!(f)?;
        }

        write!(f, "]")
    }
}

auto_ops::impl_op_ex!(+=|a: &mut Matrix, b: &Matrix| {
    detail::check_size(a.row_size(), b.row_size());
    detail::check_size(a.col_size(), b.col_size());

    for (ar, br) in a.rows.iter_mut().zip(&b.rows) {
        *ar += br;
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

    for (ar, br) in a.rows.iter_mut().zip(&b.rows) {
        *ar -= br;
    }
});

auto_ops::impl_op_ex!(-|a: &Matrix, b: &Matrix| -> Matrix {
    let mut a = a.clone();
    a -= b;
    a
});

auto_ops::impl_op_ex!(*=|a: &mut Matrix, b: Fraction| {
    for row in &mut a.rows {
        *row *= b;
    }
});

auto_ops::impl_op_ex_commutative!(*|a: Matrix, b: Fraction| -> Matrix {
    let mut a = a;
    a *= b;
    a
});

auto_ops::impl_op_ex!(*=|a: &mut Matrix, b: i32| {
    for row in &mut a.rows {
        *row *= b;
    }
});

auto_ops::impl_op_ex_commutative!(*|a: Matrix, b: i32| -> Matrix {
    let mut a = a;
    a *= b;
    a
});

auto_ops::impl_op_ex!(/=|a: &mut Matrix, b: Fraction| {
    for row in &mut a.rows {
        *row /= b;
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
                *elem = -*elem;
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
