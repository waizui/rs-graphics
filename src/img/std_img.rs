use image::Rgb;
use crate::img::img_core::*;

pub struct RgbImage<T> {
    shape: (usize, usize),
    pixels: Vec<Rgb<T>>,
}

impl<T> Image2D<Rgb<T>> for RgbImage<T> {
    fn shape(&self) -> (usize, usize) {
        self.shape
    }

    fn pixels(&self) -> &[Rgb<T>] {
        &self.pixels
    }

    fn new(shape: (usize, usize), pixels: Vec<Rgb<T>>) -> Option<Self> {
        if shape.0 * shape.1 != pixels.len() {
            None
        } else {
            Some(RgbImage { shape, pixels })
        }
    }
}
