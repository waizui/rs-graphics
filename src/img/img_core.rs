pub trait Image2D<P>: Sized {
    fn shape(&self) -> (usize, usize);
    fn pixels(&self) -> &[P];
    fn new(shape: (usize, usize), pixels: Vec<P>) -> Option<Self>;
}
