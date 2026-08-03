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
    ///
    /// An empty vector is considered a zero vector.
    pub fn is_zero(&self) -> bool {
        self.count_leading_zeros() == self.size()
    }

    /// Determine whether two vectors are orthogonal.
    ///
    /// A zero vector is orthogonal to any vector, and an empty vector is
    /// orthogonal to an empty vector.
    pub fn is_orthogonal(&self, that: &Self) -> bool {
        detail::check_size(self.size(), that.size());

        (self * that) == 0.into()
    }

    /// Determine whether two vectors are parallel.
    ///
    /// A zero vector is parallel to any vector, and an empty vector is
    /// parallel to an empty vector.
    pub fn is_parallel(&self, that: &Self) -> bool {
        detail::check_size(self.size(), that.size());

        // zero vector (incl. empty) parallel to any vector
        if self.is_zero() || that.is_zero() {
            return true;
        }

        // find the first non-zero element
        let i = self.count_leading_zeros();
        let scale = that[i] / self[i];
        // compare element-by-element
        self.elements.iter().zip(&that.elements).all(|(a, b)| *a * scale == *b)
    }

    /// Calculate the squared Euclidean norm of the vector.
    ///
    /// Unlike [`norm`](Self::norm), this uses exact rational arithmetic.
    pub fn norm_squared(&self) -> Fraction {
        self.elements.iter().fold(Fraction::new(), |acc, x| acc + *x * *x)
    }

    /// Calculate the Euclidean norm of the vector.
    ///
    /// An empty vector has norm `0.0`.
    ///
    /// The result is an `f64` because the norm may be irrational; for an
    /// exact result use [`norm_squared`](Self::norm_squared).
    pub fn norm(&self) -> f64 {
        f64::from(self.norm_squared()).sqrt()
    }

    /// Calculate the number of leading zeros of this vector.
    pub fn count_leading_zeros(&self) -> usize {
        self.elements.iter().position(|x| *x != 0.into()).unwrap_or(self.size())
    }

    /// Return the cross product of two vectors.
    ///
    /// 2D vectors are embedded in the `xy`-plane (zero `z`-component), so the
    /// result is always a 3-element vector; the 2D cross product
    /// `a[0] * b[1] - a[1] * b[0]` is the `z`-component of the result.
    ///
    /// # Panics
    /// Panics unless both vectors have size 2 or both have size 3.
    pub fn cross(a: &Self, b: &Self) -> Self {
        match (a.size(), b.size()) {
            (2, 2) => Self::from([0.into(), 0.into(), a[0] * b[1] - a[1] * b[0]]),
            (3, 3) => Self::from([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]),
            _ => panic!("Error: Incompatible dimensions for cross product."),
        }
    }

    /// Return an iterator over the elements.
    pub fn iter(&self) -> std::slice::Iter<'_, Fraction> {
        self.elements.iter()
    }

    /// Return a mutable iterator over the elements.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Fraction> {
        self.elements.iter_mut()
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

impl FromIterator<Fraction> for Vector {
    fn from_iter<T: IntoIterator<Item = Fraction>>(iter: T) -> Self {
        Self {
            elements: iter.into_iter().collect(),
        }
    }
}

impl FromIterator<i32> for Vector {
    fn from_iter<T: IntoIterator<Item = i32>>(iter: T) -> Self {
        Self {
            elements: iter.into_iter().map(Fraction::from).collect(),
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
        use std::fmt::Write;

        write!(f, "[")?;

        // calc the max width of element
        let mut buf = String::new();
        let mut width = 0;
        for i in 0..self.size() {
            write!(buf, "{}", self[i]).unwrap();
            width = width.max(buf.len());
            buf.clear();
        }

        // align right, fill with space
        for i in 0..self.size() {
            if i != 0 {
                write!(f, " ")?;
            }
            write!(buf, "{}", self[i]).unwrap();
            write!(f, "{:>width$}", buf)?;
            buf.clear();
        }

        write!(f, "]")
    }
}

auto_ops::impl_op_ex!(+=|a: &mut Vector, b: &Vector| {
    detail::check_size(a.size(), b.size());
    for (a_i, b_i) in a.elements.iter_mut().zip(&b.elements) {
        *a_i += *b_i;
    }
});

auto_ops::impl_op_ex!(+|a: &Vector, b: &Vector| -> Vector {
    let mut a = a.clone();
    a += b;
    a
});

auto_ops::impl_op_ex!(-=|a: &mut Vector, b: &Vector| {
    detail::check_size(a.size(), b.size());
    for (a_i, b_i) in a.elements.iter_mut().zip(&b.elements) {
        *a_i -= *b_i;
    }
});

auto_ops::impl_op_ex!(-|a: &Vector, b: &Vector| -> Vector {
    let mut a = a.clone();
    a -= b;
    a
});

auto_ops::impl_op_ex!(*=|a: &mut Vector, b: Fraction| {
    for x in &mut a.elements {
        *x *= b;
    }
});

auto_ops::impl_op_ex_commutative!(*|a: Vector, b: Fraction| -> Vector {
    let mut a = a;
    a *= b;
    a
});

auto_ops::impl_op_ex!(*=|a: &mut Vector, b: i32| {
    for x in &mut a.elements {
        *x *= Fraction::from(b);
    }
});

auto_ops::impl_op_ex_commutative!(*|a: Vector, b: i32| -> Vector {
    let mut a = a;
    a *= b;
    a
});

auto_ops::impl_op_ex!(/=|a: &mut Vector, b: Fraction| {
    for x in &mut a.elements {
        *x /= b;
    }
});

auto_ops::impl_op_ex!(/|a: Vector, b: Fraction| -> Vector {
    let mut a = a;
    a /= b;
    a
});

auto_ops::impl_op_ex!(*|a: &Vector, b: &Vector| -> Fraction {
    detail::check_size(a.size(), b.size());

    // empty vectors: empty sum, i.e. the dot product is 0
    let mut result = Fraction::new();
    for i in 0..a.size() {
        result += a[i] * b[i];
    }
    result
});

impl std::ops::Neg for Vector {
    type Output = Vector;

    fn neg(mut self) -> Vector {
        for elem in &mut self.elements {
            *elem = -*elem;
        }
        self
    }
}

impl std::ops::Neg for &Vector {
    type Output = Vector;

    fn neg(self) -> Vector {
        -(self.clone())
    }
}

impl IntoIterator for Vector {
    type Item = Fraction;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.elements.into_iter()
    }
}

impl<'a> IntoIterator for &'a Vector {
    type Item = &'a Fraction;
    type IntoIter = std::slice::Iter<'a, Fraction>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a> IntoIterator for &'a mut Vector {
    type Item = &'a mut Fraction;
    type IntoIter = std::slice::IterMut<'a, Fraction>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}
