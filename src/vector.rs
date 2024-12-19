use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use crate::detail;

use pyinrs::Fraction;

/// Vector with fractions as elements.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Vector {
    pub(crate) elements: Vec<Fraction>,
}

impl Vector {
    /// Create a new vector object.
    pub fn new() -> Self {
        Self { elements: Vec::new() }
    }

    /// Create an n-dimensional vector with all identical elements.
    pub fn create(n: usize, value: Fraction) -> Self {
        Self { elements: [value].repeat(n) }
    }

    /// Create an n-dimensional vector with all 0 elements.
    pub fn zeros(n: usize) -> Self {
        Self::create(n, 0.into())
    }

    /// Create an n-dimensional vector with all 1 elements.
    pub fn ones(n: usize) -> Self {
        Self::create(n, 1.into())
    }

    /// Returns the number of elements in the vector.
    pub fn size(&self) -> usize {
        self.elements.len()
    }

    /// Returns `true` if the vector contains no elements.
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Determine if it is a zero vector.
    pub fn is_zero(&self) -> bool {
        self.count_leading_zeros() == self.size()
    }

    /// Determine whether two vectors are orthogonal.
    pub fn is_orthogonal(&self, that: &Self) -> bool {
        detail::check_empty(self.size());
        detail::check_size(self.size(), that.size());

        (self * that) == 0.into()
    }

    /// Determine whether two vectors are paralle.
    pub fn is_parallel(&self, that: &Self) -> bool {
        detail::check_empty(self.size());
        detail::check_size(self.size(), that.size());

        // zero vector parallel to any vector
        if self.is_zero() || that.is_zero() {
            return true;
        }

        // find the first non-zero element
        let i = self.count_leading_zeros();
        // calc the scale factor
        let scale = that[i] / self[i];
        // if equal after scale-up, then parallel
        self * scale == *that
    }

    /// Calculate the norm (abs) of the vector.
    pub fn norm(&self) -> f64 {
        detail::check_empty(self.size());

        let mut norm = 0.0;
        for i in 0..self.size() {
            norm += f64::from(self.elements[i] * self.elements[i]);
        }
        norm.sqrt()
    }

    /// Calculate the number of leading zeros of this vector.
    pub fn count_leading_zeros(&self) -> usize {
        detail::check_empty(self.size());

        let mut lz: usize = 0;
        while self.elements[lz] == 0.into() {
            lz += 1;
            if lz == self.size() {
                break;
            }
        }
        lz
    }

    /// Return the cross product of two vectors.
    pub fn cross(a: &Self, b: &Self) -> Self {
        if a.size() == 2 && b.size() == 2 {
            Self::from([a[0] * b[1] - a[1] * b[0]])
        } else if a.size() == 3 && b.size() == 3 {
            return Self::from([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]);
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
        Self {
            elements: Vec::from(value.map(Fraction::from)),
        }
    }
}

impl From<Vec<Fraction>> for Vector {
    fn from(value: Vec<Fraction>) -> Self {
        Self { elements: value }
    }
}

impl From<Vec<i32>> for Vector {
    fn from(value: Vec<i32>) -> Self {
        Self {
            elements: value.into_iter().map(Fraction::from).collect(),
        }
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
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;

        // calc the max width of element
        let mut width = 0;
        for i in 0..self.size() {
            width = width.max(format!("{}", self[i]).len());
        }

        // align right, fill with space
        for i in 0..self.size() {
            if i != 0 {
                write!(f, " ")?;
            }
            write!(f, "{:>width$}", format!("{}", self[i]))?;
        }

        write!(f, "]")
    }
}

auto_ops::impl_op_ex!(+=|a: &mut Vector, b: &Vector| {
    detail::check_size(a.size(), b.size());
    for i in 0..a.size() {
        a[i] += b[i];
    };
});

auto_ops::impl_op_ex!(+|a: &Vector, b: &Vector| -> Vector {
    let mut a = a.clone();
    a += b;
    a
});

auto_ops::impl_op_ex!(-=|a: &mut Vector, b: &Vector| {
    detail::check_size(a.size(), b.size());
    for i in 0..a.size() {
        a[i] -= b[i];
    };
});

auto_ops::impl_op_ex!(-|a: &Vector, b: &Vector| -> Vector {
    let mut a = a.clone();
    a -= b;
    a
});

auto_ops::impl_op_ex!(*=|a: &mut Vector, b: Fraction| {
    for i in 0..a.size() {
        a[i] *= b;
    }
});

auto_ops::impl_op_ex_commutative!(*|a: &Vector, b: Fraction| -> Vector {
    let mut a = a.clone();
    a *= b;
    a
});

auto_ops::impl_op_ex!(*=|a: &mut Vector, b: i32| {
    for i in 0..a.size() {
        a[i] *= Fraction::from(b);
    }
});

auto_ops::impl_op_ex_commutative!(*|a: &Vector, b: i32| -> Vector {
    let mut a = a.clone();
    a *= b;
    a
});

auto_ops::impl_op_ex!(*|a: &Vector, b: &Vector| -> Fraction {
    detail::check_empty(a.size());
    detail::check_size(a.size(), b.size());

    let mut result = Fraction::new();
    for i in 0..a.size() {
        result += a[i] * b[i];
    }
    result
});

impl IntoIterator for Vector {
    type Item = Fraction;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.elements.into_iter()
    }
}
