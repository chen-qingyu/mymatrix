use std::{
    fmt::Display,
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
};

use crate::utility;

use pyinrs::Fraction;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Vector {
    pub(crate) elements: Vec<Fraction>,
}

impl Vector {
    /// Construct a new vector object.
    pub fn new() -> Vector {
        Self { elements: Vec::new() }
    }

    /// Construct a vector with n identical elements.
    pub fn create(n: usize, value: Fraction) -> Vector {
        Self { elements: [value].repeat(n) }
    }

    /// Returns the number of elements in the vector.
    pub fn size(&self) -> usize {
        self.elements.len()
    }

    /// Returns `true` if the vector contains no elements.
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Returns an iterator over the vector.
    pub fn iter(&self) -> std::slice::Iter<Fraction> {
        self.elements.iter()
    }

    /// Calculate the length of the vector.
    pub fn length(&self) -> f64 {
        utility::check_empty(self.size());

        let mut length = 0.0;
        for i in 0..self.size() {
            length += f64::from(self.elements[i] * self.elements[i]);
        }
        length.sqrt()
    }

    /// Calculate the number of leading zeros for this vector.
    pub fn count_leading_zeros(&self) -> usize {
        utility::check_empty(self.size());

        let mut lz: usize = 0;
        while self.elements[lz] == 0.into() {
            lz += 1;
            if lz == self.size() {
                break;
            }
        }
        lz
    }

    /// Determine if it is a zero vector.
    pub fn is_zero(&self) -> bool {
        self.count_leading_zeros() == self.size()
    }

    /// Determine whether two vectors are orthogonal.
    pub fn is_orthogonal(&self, vector: &Self) -> bool {
        utility::check_empty(self.size());
        utility::check_size(self.size(), vector.size());

        (self * vector) == 0.into()
    }

    /// Determine whether two vectors are paralle.
    pub fn is_parallel(&self, vector: &Self) -> bool {
        utility::check_empty(self.size());
        utility::check_size(self.size(), vector.size());

        f64::from(self * vector).abs() == self.length() * vector.length()
    }

    /// Return the cross product of two vectors.
    pub fn cross(a: &Self, b: &Self) -> Self {
        if a.size() == 2 && b.size() == 2 {
            Vector::from([a[0] * b[1] - a[1] * b[0]])
        } else if a.size() == 3 && b.size() == 3 {
            return Vector::from([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]);
        } else {
            panic!("Error: Incompatible dimensions for cross product.");
        }
    }
}

impl<const N: usize> From<[Fraction; N]> for Vector {
    fn from(value: [Fraction; N]) -> Self {
        Self { elements: Vec::from(value) }
    }
}

impl<const N: usize> From<[i32; N]> for Vector {
    fn from(value: [i32; N]) -> Self {
        let v: [Fraction; N] = value.map(|x| x.into());
        Self { elements: Vec::from(v) }
    }
}

impl Index<usize> for Vector {
    type Output = Fraction;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elements[index]
    }
}

impl IndexMut<usize> for Vector {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elements[index]
    }
}

impl Display for Vector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for i in 0..self.elements.len() {
            if i > 0 {
                write!(f, " ")?;
            }
            write!(f, "{}", self.elements[i])?;
        }
        write!(f, "]")
    }
}

impl Add<&Vector> for &Vector {
    type Output = Vector;

    fn add(self, rhs: &Vector) -> Self::Output {
        utility::check_size(self.size(), rhs.size());

        let mut result = self.clone();
        for i in 0..self.size() {
            result[i] += rhs[i];
        }
        result
    }
}

impl Sub<&Vector> for &Vector {
    type Output = Vector;

    fn sub(self, rhs: &Vector) -> Self::Output {
        utility::check_size(self.size(), rhs.size());

        let mut result = self.clone();
        for i in 0..self.size() {
            result[i] -= rhs[i];
        }
        result
    }
}

impl Mul<&Vector> for &Vector {
    type Output = Fraction;

    fn mul(self, rhs: &Vector) -> Self::Output {
        utility::check_empty(self.size());
        utility::check_size(self.size(), rhs.size());

        let mut result = 0.into();
        for i in 0..self.size() {
            result += self[i] * rhs[i];
        }
        result
    }
}

impl Mul<Fraction> for &Vector {
    type Output = Vector;

    fn mul(self, rhs: Fraction) -> Self::Output {
        let mut result = self.clone();
        for i in 0..self.size() {
            result[i] *= rhs;
        }
        result
    }
}

impl Add for Vector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl Sub for Vector {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl Mul for Vector {
    type Output = Fraction;

    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}

impl Mul<Fraction> for Vector {
    type Output = Self;

    fn mul(mut self, rhs: Fraction) -> Self::Output {
        for i in 0..self.size() {
            self[i] *= rhs;
        }
        self
    }
}

impl AddAssign<&Vector> for Vector {
    fn add_assign(&mut self, rhs: &Vector) {
        *self = &*self + rhs;
    }
}

impl SubAssign<&Vector> for Vector {
    fn sub_assign(&mut self, rhs: &Vector) {
        *self = &*self - rhs;
    }
}

impl MulAssign<Fraction> for Vector {
    fn mul_assign(&mut self, rhs: Fraction) {
        *self = &*self * rhs;
    }
}
