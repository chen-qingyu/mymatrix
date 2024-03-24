#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Vector<T> {
    elements: Vec<T>,
}

impl<T: Clone> Vector<T> {
    /// Construct a new vector object.
    pub fn new() -> Self {
        Self { elements: [].to_vec() }
    }
}
