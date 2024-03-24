use std::{
    fmt::Display,
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
    usize,
};

use crate::{utility, Vector};

use pyinrs::Fraction;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Hash, Default)]
pub struct Matrix {
    rows: Vec<Vector>,
}

impl Matrix {
    /// Construct a new matrix object.
    pub fn new() -> Matrix {
        Matrix { rows: Vec::new() }
    }

    /// Construct a matrix with row x col identical elements.
    pub fn create(row: usize, col: usize, value: Fraction) -> Matrix {
        let mut m = Matrix::new();
        for _ in 0..row {
            m.rows.push(Vector::create(col, value));
        }
        m
    }

    /// Generate an n-order unit matrix.
    pub fn eye(n: usize) -> Matrix {
        let mut m = Matrix::create(n, n, 0.into());
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

    /// Calculate the rank of this matrix.
    pub fn rank(&self) -> usize {
        todo!()
    }

    /// Compute the determinant of this matrix.
    pub fn det(&self) -> Fraction {
        todo!()
    }

    /// Compute the inverse of this matrix.
    pub fn inv(&self) -> Matrix {
        todo!()
    }

    /// Expand this matrix by rows.
    pub fn append_row(&mut self, mut matrix: Matrix) -> &Matrix {
        utility::check_size(self.col_size(), matrix.col_size());

        self.rows.append(&mut matrix.rows);
        self
    }

    /// Expand this matrix by columns.
    pub fn append_col(&mut self, mut matrix: Matrix) -> &Matrix {
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
        for (e, &value) in self.clone()[j].iter().enumerate() {
            self[i][e] += value * k;
        }
        self
    }

    /// Transform this matrix to general row echelon form.
    pub fn transform_row_echelon(&mut self) -> &Self {
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

        // step 2: Transform to the row echelon form. It's so elegant, I'm a genius haha.
        self.rows.sort_by_key(|r| r.count_leading_zeros());
        self
    }

    /// Split this matrix by rows.
    pub fn split_row(&self, n: usize) -> (Matrix, Matrix) {
        utility::check_bounds(n, 0, self.row_size());

        let mut first = Matrix::new();
        let mut second = Matrix::new();
        first.rows = self.rows[0..n].to_vec();
        second.rows = self.rows[n..].to_vec();

        return (first, second);
    }

    /// Split this matrix by columns.
    pub fn split_col(&self, n: usize) -> (Matrix, Matrix) {
        utility::check_bounds(n, 0, self.col_size());

        let mut first = Matrix::new();
        let mut second = Matrix::new();
        first.rows.reserve(self.row_size());
        second.rows.reserve(self.row_size());
        for r in 0..self.row_size() {
            first.rows[r].elements = self.rows[r].elements[0..n].to_vec();
            second.rows[r].elements = self.rows[r].elements[n..].to_vec();
        }

        return (first, second);
    }

    /// Returns the transpose of the matrix.
    pub fn transpose(&self) -> Matrix {
        let mut result = Matrix::create(self.col_size(), self.row_size(), 0.into());

        for i in 0..self.row_size() {
            for j in 0..self.col_size() {
                result[j][i] = self[i][j];
            }
        }
        return result;
    }

    /// Return the product of two matrices.
    pub fn dot(a: &Matrix, b: &Matrix) -> Matrix {
        utility::check_size(a.col_size(), b.row_size());

        let mut result = Matrix::create(a.row_size(), b.col_size(), 0.into());
        let mt = b.transpose();
        for r in 0..a.row_size() {
            for c in 0..b.col_size() {
                result[r][c] = Vector::dot(&a[r], &mt[c]);
            }
        }
        return result;
    }
}

impl<const N: usize> From<[Vector; N]> for Matrix {
    fn from(value: [Vector; N]) -> Self {
        Self { rows: Vec::from(value) }
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
        for i in 0..self.row_size() {
            write!(f, "{}\n", self.rows[i])?;
        }
        Ok(())
    }
}

impl AddAssign for Matrix {
    fn add_assign(&mut self, rhs: Self) {
        utility::check_size(self.row_size(), rhs.row_size());
        utility::check_size(self.col_size(), rhs.col_size());

        for r in 0..self.row_size() {
            self.rows[r] += &rhs[r];
        }
    }
}

impl SubAssign for Matrix {
    fn sub_assign(&mut self, rhs: Self) {
        utility::check_size(self.row_size(), rhs.row_size());
        utility::check_size(self.col_size(), rhs.col_size());

        for r in 0..self.row_size() {
            self.rows[r] -= &rhs[r];
        }
    }
}

impl MulAssign for Matrix {
    fn mul_assign(&mut self, rhs: Self) {
        utility::check_size(self.row_size(), rhs.row_size());
        utility::check_size(self.col_size(), rhs.col_size());

        for r in 0..self.row_size() {
            self.rows[r] *= &rhs[r];
        }
    }
}

impl MulAssign<Fraction> for Matrix {
    fn mul_assign(&mut self, rhs: Fraction) {
        for r in 0..self.row_size() {
            self.rows[r] *= rhs;
        }
    }
}

impl Add for Matrix {
    type Output = Matrix;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl Sub for Matrix {
    type Output = Matrix;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= rhs;
        self
    }
}

impl Mul for Matrix {
    type Output = Matrix;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl Mul<Fraction> for Matrix {
    type Output = Matrix;

    fn mul(mut self, rhs: Fraction) -> Self::Output {
        self *= rhs;
        self
    }
}
