# MyMatrix

_My simple matrix library that can perform fraction operations._

## 1. Attribute

- Name: MyMatrix
- Goal: Provide a simple matrix library that can perform fraction operations
- Module: Fraction, Vector, Matrix
- Test: Using [rstest](https://crates.io/crates/rstest) for unit tests and ensure all tests passed
- Security: There is no `unsafe` code block

## 2. Usage

To use it, add the following lines to your `Cargo.toml` file:

```toml
[dependencies]
mymatrix = "1.1"
```

Some simple examples:

```rust
use mymatrix::{Fraction, Vector, Matrix};

// Vector dot product
Vector::from([1, 2, 3]) * Vector::from([4, 5, 6]); // 32
// Vector cross product
Vector::cross(&[1, 2, 3].into(), &[4, 5, 6].into()); // [-3  6 -3]
// Vector scalar product
Vector::from([1, 2, 3]) * Fraction::from((2, 5)); // [2/5 4/5 6/5]

// Matrix rank
Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).rank(); // 3
// Matrix determinant
Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).det(); // 27
// Matrix inversion
Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).inv().unwrap();
/*
[
-16/9   8/9  -1/9
 14/9  -7/9   2/9
 -1/9   2/9  -1/9
]
*/

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

let A = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]);
assert_eq!(A.adj(), A.det() * A.inv().unwrap()); //  A.adj  = |A| * A.inv
assert_eq!(A.adj().det(), A.det() * A.det());    // |A.adj| = |A|^(n-1)

// Matrix * Vector
Matrix::from([[1, 2], [3, 4]]) * Vector::from([1, 2]); // [5, 11]

// Solve linear system Ax = b
let a = Matrix::from([[2, 3], [4, 5]]);
let b = Vector::from([7, 13]);
let x = a.solve(&b).unwrap(); // [2, 1]

// General solution: (particular solution, null space basis)
let m = Matrix::from([[1, 1], [1, 1]]);
let (xp, null) = m.general_solution(&Vector::from([2, 2])).unwrap();
// xp = [2, 0], null = [-1, 1]  =>  every solution is xp + c * null[0]

// Matrix power
Matrix::from([[1, 2], [3, 4]]).pow(3);
/*
[
37  54
81 118
]
*/

// Subspace bases
let m = Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
m.row_space();  // row space basis:    [1 2 3; 0 -3 -6]
m.col_space();  // column space basis: [1 4 7; 2 5 8]
m.null_space(); // null space basis:   [1 -2 1]

// Characteristic polynomial: det(λI - A) = λ³ - 9λ² + 26λ - 24
Matrix::from([[2, 0, 0], [0, 3, 0], [0, 0, 4]]).characteristic_polynomial(); // [1 -9 26 -24]

// Moore-Penrose pseudo-inverse (works for singular / non-square matrices)
Matrix::from([[1, 2], [2, 4]]).pseudo_inverse();
/*
[
1/25 2/25
2/25 4/25
]
*/

// Column accessor
Matrix::from([[1, 2, 3], [4, 5, 6]]).col(1); // [2 5]

// Kronecker product and Hadamard (element-wise) product
Matrix::from([[1, 2], [3, 4]]).kron(&Matrix::identity(2));
Matrix::from([[1, 2], [3, 4]]).hadamard(&Matrix::ones(2, 2)); // [1 2; 3 4]

// Diagonal extraction and construction
Matrix::identity(3).diag(); // [1 1 1]
Matrix::from_diagonal(&Vector::from([2, 3, 4]));
```
