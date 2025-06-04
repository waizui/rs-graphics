pub mod operation;

/// row major bitmap
pub struct Image<T> {
    pub shape: (usize, usize),
    pub data: Vec<T>,
}
