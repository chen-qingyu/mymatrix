# MyMatrix

_My simple matrix library that can perform fraction operations._

## 1. Attribute

- Name: MyMatrix.
- Language: Rust, requires version rustc >= `1.75.0`.
- Goal: Write a simple matrix library that can perform fraction operations.
- Module: Vector, Matrix
- Style: Follow Rust's official recommended style.
- Test: Using [rstest](https://crates.io/crates/rstest) for unit testing and ensure that all tests passed.
- Security: There is no `unsafe` code block.
- Document: Using `cargo doc --open` to open documents.

## 2. Usage

To use it, add the following lines to your `Cargo.toml` file:

```toml
[dependencies]
mymatrix = "0"
```

Some simple examples:

```rust
use mymatrix::*;

// Vector dot product
print!("{}", Vector::from([1, 2, 3]) * Vector::from([4, 5, 6])); // 32
// Vector cross product
print!("{}", Vector::cross(&[1, 2, 3].into(), &[4, 5, 6].into())); // [-3 6 -3]
// Vector scalar product
print!("{}", Vector::from([1, 2]) * pyinrs::Fraction::from((2, 5))); // [2/5 4/5]

// Matrix rank
print!("{}", Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).rank()); // 2
// Matrix determinant
print!("{}", Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).det().unwrap()); // 27
// Matrix inversion
print!("{}", Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 0]]).inv().unwrap());
/*
[[
-16/9 8/9 -1/9;
14/9 -7/9 2/9;
-1/9 2/9 -1/9;
]]
*/

let a = Matrix::from([[1, 2], [3, 4]]);
let b = Matrix::create(2, 2, 0.into());
let c = Matrix::create(2, 2, 1.into());
let d = Matrix::eye(2);

print!("{}", ((a + b) * (c + d)).inv().unwrap());
/*
[[
-11/6 5/6;
5/3 -2/3;
]]
*/
```
