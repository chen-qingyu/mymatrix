#![doc = include_str!("../readme.md")]

mod utility;

mod matrix;
mod vector;

pub use matrix::Matrix;
pub use pyinrs::Fraction;
pub use vector::Vector;
