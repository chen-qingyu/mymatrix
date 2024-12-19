#![doc = include_str!("../readme.md")]

mod detail;

mod matrix;
mod vector;

pub use matrix::Matrix;
pub use pyinrs::Fraction;
pub use vector::Vector;
