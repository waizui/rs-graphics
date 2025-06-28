use rayon::prelude::*;

use crate::img::img_core::*;

/// image operations
pub trait ImageOp<P>: Image2D<P> + Sync
where
    P: Copy + Send + Sync,
{
    /// stitch two images horizontally
    fn stitch_hor(&self, other: &Self, defval: P) -> Option<Self> {
        let (w0, h0) = self.shape();
        let (w1, h1) = other.shape();
        let (w, h) = (w0 + w1, h0.max(h1));

        let mut pixels = vec![defval; w * h];

        pixels.par_iter_mut().enumerate().for_each(|(ipix, pix)| {
            let x = ipix % w;
            let y = ipix / w;

            if x < w0 && y < h0 {
                // pixel belongs to the first image
                *pix = self.pixels()[y * w0 + x];
            } else if x >= w0 && y < h1 {
                let x1 = x - w0;
                *pix = other.pixels()[y * w1 + x1];
            } else {
                *pix = defval;
            }
        });

        Image2D::new((w, h), pixels)
    }

    /// stitch two images vertically
    fn stitch_ver(&self, other: &Self, defval: P) -> Option<Self> {
        if self.shape() == other.shape() {
            return Image2D::new(
                (self.shape().0, other.shape().1 * 2),
                [self.pixels(), other.pixels()].concat(),
            );
        }

        let (w0, h0) = self.shape();
        let (w1, h1) = other.shape();
        let (w, h) = (w0.max(w1), h0 + h1);

        let mut pixels = vec![defval; w * h];
        pixels.par_iter_mut().enumerate().for_each(|(ipix, pix)| {
            let x = ipix % w;
            let y = ipix / w;

            if x < w0 && y < h0 {
                *pix = self.pixels()[y * w0 + x];
            } else if x < w1 && y >= h0 && y < h {
                let y1 = y - h0;
                *pix = other.pixels()[y1 * w1 + x];
            } else {
                *pix = defval;
            }
        });

        Image2D::new((w, h), pixels)
    }

    fn stitch_hor_mult(imgs: &[Self], defval: P) -> Option<Self> {
        let mut iter = imgs.iter();
        let n = iter.next()?;
        let first = Self::new(n.shape(), n.pixels().to_vec())?;
        iter.try_fold(first, |acc, img| acc.stitch_hor(img, defval))
    }
}

/// extent Image2D abilities
impl<P, T> ImageOp<P> for T
where
    P: Copy + Send + Sync,
    T: Image2D<P> + Sync,
{
}
