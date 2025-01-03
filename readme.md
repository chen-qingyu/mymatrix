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
mymatrix = "0"
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
```
