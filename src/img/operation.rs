use super::Image;
use rayon::prelude::*;


/// stitch two images horizontally
pub fn stitch_hor<T>(img0: &Image<T>, img1: &Image<T>, defval: T) -> Image<T>
where
    T: Copy + Send + Sync,
{
    let (w, h) = {
        let (w0, h0) = img0.shape;
        let (w1, h1) = img1.shape;
        (w0 + w1, h0.max(h1))
    };

    let task = |ipix: usize, pix: &mut T| {
        let (w0, h0) = img0.shape;
        let (w1, h1) = img1.shape;

        let x = ipix % w;
        let y = ipix / w;

        if x < w0 && y < h0 {
            // pixel belongs to the first image
            *pix = img0.data[y * w0 + x];
        } else if x >= w0 && y < h1 {
            let x1 = x - w0;
            *pix = img1.data[y * w1 + x1];
        } else {
            *pix = defval;
        }
    };

    let mut data = vec![defval; w * h];
    data.par_iter_mut()
        .enumerate()
        .for_each(|(ipix, pix)| task(ipix, pix));

    Image {
        shape: (w, h),
        data,
    }
}

pub fn stitch_hor_mult<T>(imgs: &[Image<T>], defval: T) -> Image<T>
where
    T: Copy + Send + Sync,
{
    let mut res = Image {
        shape: (1, 1),
        data: vec![defval],
    };

    for img in imgs.iter() {
        res = stitch_hor(&res, img, defval);
    }

    res
}

/// stitch two images vertically
pub fn stitch_ver<T>(img0: &Image<T>, img1: &Image<T>, defval: T) -> Image<T>
where
    T: Copy + Send + Sync,
{
    if img0.shape == img1.shape {
        Image {
            shape: (img0.shape.0, img0.shape.1 * 2),
            data: [img0.data.clone(), img1.data.clone()].concat(),
        }
    } else {
        let (w, h) = {
            let (w0, h0) = img0.shape;
            let (w1, h1) = img1.shape;
            (w0.max(w1), h0 + h1)
        };

        let task = |ipix: usize, pix: &mut T| {
            let (w0, h0) = img0.shape;
            let (w1, _h1) = img1.shape;

            let x = ipix % w;
            let y = ipix / w;

            if x < w0 && y < h0 {
                // pixel belongs to the first image
                *pix = img0.data[y * w0 + x];
            } else if x < w1 && y >= h0 && y < h {
                let y1 = y - h0;
                *pix = img1.data[y1 * w1 + x];
            } else {
                *pix = defval;
            }
        };

        let mut data = vec![defval; w * h];
        data.par_iter_mut()
            .enumerate()
            .for_each(|(ipix, pix)| task(ipix, pix));

        Image {
            shape: (w, h),
            data,
        }
    }
}
