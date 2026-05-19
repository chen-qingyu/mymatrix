#![doc = include_str!("../readme.md")]

mod detail;

mod error;
mod matrix;
mod vector;

pub use error::MatrixError;
pub use matrix::Matrix;
pub use pyinrs::Fraction;
pub use vector::Vector;
