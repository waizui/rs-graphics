use crate::img::img_core::*;

pub type RgbPixel<T> = image::Rgb<T>;

pub struct RgbImage<T> {
    shape: (usize, usize),
    pixels: Vec<RgbPixel<T>>,
}

impl<T> Image2D<RgbPixel<T>> for RgbImage<T> {
    fn shape(&self) -> (usize, usize) {
        self.shape
    }

    fn pixels(&self) -> &[RgbPixel<T>] {
        &self.pixels
    }

    fn new(shape: (usize, usize), pixels: Vec<RgbPixel<T>>) -> Option<Self> {
        if shape.0 * shape.1 != pixels.len() {
            None
        } else {
            Some(RgbImage { shape, pixels })
        }
    }
}
