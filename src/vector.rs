use std::{
    fmt::Display,
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
};

use crate::utility;

use pyinrs::Fraction;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Hash, Default)]
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

        Self::dot(self, vector) == 0.into()
    }

    /// Determine whether two vectors are paralle.
    pub fn is_parallel(&self, vector: &Self) -> bool {
        utility::check_empty(self.size());
        utility::check_size(self.size(), vector.size());

        f64::from(Self::dot(self, vector)).abs() == self.length() * vector.length()
    }
    /// Return the dot product (scalar product, inner product) of two vectors.
    pub fn dot(a: &Self, b: &Self) -> Fraction {
        utility::check_empty(a.size());
        utility::check_size(a.size(), b.size());

        let mut result = 0.into();
        for i in 0..a.size() {
            result += a[i] * b[i];
        }
        result
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

impl AddAssign<&Vector> for Vector {
    fn add_assign(&mut self, rhs: &Vector) {
        utility::check_empty(self.size());
        utility::check_size(self.size(), rhs.size());

        for i in 0..self.size() {
            self[i] += rhs[i];
        }
    }
}

impl SubAssign<&Vector> for Vector {
    fn sub_assign(&mut self, rhs: &Vector) {
        utility::check_empty(self.size());
        utility::check_size(self.size(), rhs.size());

        for i in 0..self.size() {
            self[i] -= rhs[i];
        }
    }
}

impl MulAssign<&Vector> for Vector {
    fn mul_assign(&mut self, rhs: &Vector) {
        utility::check_empty(self.size());
        utility::check_size(self.size(), rhs.size());

        for i in 0..self.size() {
            self[i] *= rhs[i];
        }
    }
}

impl MulAssign<Fraction> for Vector {
    fn mul_assign(&mut self, rhs: Fraction) {
        utility::check_empty(self.size());

        for i in 0..self.size() {
            self[i] *= rhs;
        }
    }
}

impl Add for Vector {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += &rhs;
        self
    }
}

impl Sub for Vector {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self -= &rhs;
        self
    }
}

impl Mul for Vector {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= &rhs;
        self
    }
}

impl Mul<Fraction> for Vector {
    type Output = Self;

    fn mul(mut self, rhs: Fraction) -> Self::Output {
        self *= rhs;
        self
    }
}
