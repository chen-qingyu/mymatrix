# MyMatrix

_My simple matrix library that can perform fraction operations._

## 1. Attribute

- **Name**: MyMatrix
- **Goal**: Provide a simple matrix library that can perform fraction operations
- **Module**: `Fraction`, `Vector`, `Matrix`
- **Exactness**: Every operation uses exact rational arithmetic (`pyinrs::Fraction`), so there is no floating-point error — the only exception is `Vector::norm()`, which returns an `f64`
- **Error handling**: structurally invalid input (non-square, mismatched dimensions, out-of-bounds) panics; mathematically invalid input (singular, inconsistent, not positive definite) is returned as `Result<_, MatrixError>`
- **Test**: Using [rstest](https://crates.io/crates/rstest) for unit tests and ensure all tests passed
- **Security**: There is no `unsafe` code block

## 2. Feature Overview

- **Vector**: dot / cross / norm / norm_squared / zero / orthogonal / parallel
- **Matrix arithmetic**: add, subtract, multiply, scalar operations, transpose, rotations, Kronecker & Hadamard products, power
- **Linear algebra**: determinant, inverse, adjugate, minor / cofactor, rank, row / column / null space bases, Moore-Penrose pseudo-inverse
- **Decompositions**: LDLᵀ (exact-rational Cholesky) and LU (Doolittle, with row pivoting)
- **Polynomials**: characteristic polynomial (Newton's identities)
- **Linear systems**: `solve` (unique solution) and `general_solution` (particular solution + null space basis)
- **Structure**: row echelon / reduced row echelon forms, row / column split & expand, elementary row operations, diagonal extraction / construction

## 3. Usage

To use it, add the following lines to your `Cargo.toml` file:

```toml
[dependencies]
mymatrix = "1.2"
```

### 3.1 Vector operations

```rust
use mymatrix::{Fraction, Vector};

// dot product (returns a Fraction)
Vector::from([1, 2, 3]) * Vector::from([4, 5, 6]); // 32

// cross product; 2D vectors are embedded in the xy-plane, so the
// result is always a 3-element vector
Vector::cross(&[1, 2].into(), &[3, 4].into()); // [0 0 -2]

// scalar multiplication (by Fraction or i32)
Vector::from([1, 2, 3]) * Fraction::from((2, 5)); // [2/5 4/5 6/5]

// queries
Vector::from([3, 4]).norm();                                // 5.0
Vector::from([3, 4]).norm_squared();                        // 25 (exact)
Vector::from([0, 0]).is_zero();                             // true
Vector::from([1, 1]).is_orthogonal(&Vector::from([1, -1])); // true
Vector::from([1, 1]).is_parallel(&Vector::from([2, 2]));    // true
```

### 3.2 Construction

```rust
use mymatrix::{Fraction, Matrix, Vector};

// from arrays, Vecs, or iterators
Matrix::from([[1, 2], [3, 4]]);
let _: Matrix = [vec![1, 2], vec![3, 4]].into_iter().collect();
let _: Vector = (0..3).map(Fraction::from).collect();

// constant-filled matrices
Matrix::zeros(2, 3);
Matrix::ones(2, 3);
Matrix::identity(3);
Matrix::create(2, 2, Fraction::from((1, 2)));

// diagonal construction / extraction
Matrix::from_diagonal(&Vector::from([2, 3, 4]));
Matrix::identity(3).diag(); // [1 1 1]
```

### 3.3 Matrix arithmetic

```rust
use mymatrix::{Matrix, Vector};

let a = Matrix::from([[1, 2], [3, 4]]);
let b = Matrix::zeros(2, 2);
let c = Matrix::ones(2, 2);
let d = Matrix::identity(2);

((a + b) * (c + d)).inv().unwrap();
/*
[
-11/6   5/6
  5/3  -2/3
]
*/

// matrix-vector product (column-vector convention)
Matrix::from([[1, 2], [3, 4]]) * Vector::from([1, 2]); // [5, 11]

// Kronecker and Hadamard (element-wise) products
Matrix::from([[1, 2], [3, 4]]).kron(&Matrix::identity(2));
Matrix::from([[1, 2], [3, 4]]).hadamard(&Matrix::ones(2, 2)); // [1 2; 3 4]

// power (binary exponentiation; pow(0) = identity)
Matrix::from([[1, 2], [3, 4]]).pow(3);
/*
[
37  54
81 118
]
*/
```

### 3.4 Queries and structure

```rust
use mymatrix::Matrix;

let m = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);

m.trace();         // 15
m.is_square();     // true
m.is_symmetric();  // false
m.is_upper();      // false
m.is_lower();      // false
m.is_diagonal();   // false
Matrix::zeros(2, 3).is_zero(); // true

m.transpose();
m.rotate_left();   // 90 degrees counter-clockwise
m.rotate_right();  // 90 degrees clockwise

Matrix::from([[1, 2, 3], [4, 5, 6]]).col(1); // [2 5]
```

### 3.5 Determinant, inverse, adjugate

```rust
use mymatrix::Matrix;

let m = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]);

m.rank(); // 3
m.det();  // 27
m.inv().unwrap();
/*
[
-16/9   8/9  -1/9
 14/9  -7/9   2/9
 -1/9   2/9  -1/9
]
*/

let A = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]);
assert_eq!(A.adj(), A.det() * A.inv().unwrap()); //  A.adj  = |A| * A.inv
assert_eq!(A.adj().det(), A.det() * A.det());    // |A.adj| = |A|^(n-1)
```

### 3.6 Solving linear systems

```rust
use mymatrix::{Matrix, MatrixError, Vector};

// unique solution
let a = Matrix::from([[2, 3], [4, 5]]);
let b = Vector::from([7, 13]);
let x = a.solve(&b).unwrap(); // [2, 1]

// general solution: (particular solution, null space basis)
let m = Matrix::from([[1, 1], [1, 1]]);
let (xp, null) = m.general_solution(&Vector::from([2, 2])).unwrap();
// xp = [2, 0], null = [-1, 1]  =>  every solution is xp + c * null[0]

// errors distinguish "no solution" from "no unique solution"
let inconsistent = Matrix::from([[1, 1], [1, 1]]);
assert_eq!(inconsistent.solve(&Vector::from([2, 3])), Err(MatrixError::Inconsistent));
let singular = Matrix::from([[1, 2], [2, 4]]);
assert_eq!(singular.solve(&Vector::from([3, 6])), Err(MatrixError::Singular));
```

### 3.7 Subspaces and decompositions

```rust
use mymatrix::Matrix;

// subspace bases (each basis vector is a row of the result)
let m = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
m.row_space();  // row space basis:    [1 2 3; 0 -3 -6]
m.col_space();  // column space basis: [1 4 7; 2 5 8]
m.null_space(); // null space basis:   [1 -2 1]

// LDL^T decomposition (exact-rational Cholesky): A = L * D * L^T
let (l, d) = Matrix::from([[4, 2], [2, 3]]).cholesky().unwrap();
// l = [1 0; 1/2 1], d = [4 2]

// LU decomposition (Doolittle, with row pivoting)
Matrix::from([[2, 3, 1], [4, 7, 1], [6, 7, 3]]).lu_decomposition().unwrap();

// characteristic polynomial: det(λI - A) = λ³ - 9λ² + 26λ - 24
Matrix::from([[2, 0, 0], [0, 3, 0], [0, 0, 4]]).characteristic_polynomial(); // [1 -9 26 -24]

// Moore-Penrose pseudo-inverse (also works for singular / non-square)
Matrix::from([[1, 2], [2, 4]]).pseudo_inverse();
/*
[
1/25 2/25
2/25 4/25
]
*/
```

### 3.8 Echelon forms and row operations

```rust
use mymatrix::Matrix;

let m = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
m.row_echelon_form();    // [1 2 3; 0 -3 -6; 0 0 0]
m.row_canonical_form();  // [1 0 -1; 0 1 2; 0 0 0]

// elementary row operations (chainable)
let mut m = Matrix::from([[1, 2], [3, 4]]);
m.e_row_swap(0, 1);
m.e_scalar_multiplication(1, 2.into());
m.e_row_sum(0, 1, (-1).into());

// split / expand rows or columns
let m = Matrix::from([[1, 2], [3, 4], [5, 6]]);
m.split_row(1);
m.split_col(1);
```

### 3.9 Iterating

```rust
use mymatrix::{Fraction, Matrix, Vector};

let v = Vector::from([1, 2, 3]);
let m = Matrix::from([[1, 2], [3, 4]]);

// by value, by reference, or by mutable reference
for x in &v {
    let _: &Fraction = x;
}
for row in &m {
    let _: &Vector = row;
}
let mut w = v.clone();
for x in &mut w {
    *x *= Fraction::from(2);
}
```

## 4. Notes

- **Column-vector convention**: `Matrix * Vector` treats the vector as a column; `Vector * Vector` is the dot product and returns a `Fraction`.
- **Empty vector** is the zero vector of R⁰: `is_zero()` is `true`, the empty dot product is `0`, and `norm()` is `0`.
- **Empty matrix**: `det(0x0) = 1` (empty-product convention) and the inverse of the empty matrix is itself.
