use del_raycast_core;
use rs_sampler::{envmap::pixel2texpair, haltonsampler::HaltonSampler, sampler::Sampler};

type Real = f32;
type Rgb = image::Rgb<Real>;

fn main() {
    use image::Pixel;
    use itertools::Itertools;
    use rayon::iter::IndexedParallelIterator;
    use rayon::iter::IntoParallelRefMutIterator;
    use rayon::iter::ParallelIterator;

    let w = 1024;
    let h = 1024;

    let len_rad = 0.0125;
    let focal_dis = 0.05;
    let fov = 20.0;

    let transform_cam_lcl2glbl = del_geo_core::mat4_col_major::from_translate(&[0., 0., -5.]);
    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    let iter = |i_pix: usize, pix: &mut Rgb| {
        let iw = i_pix % w;
        // top-left to bottom-left
        let ih = w - 1 - i_pix / w;
        let tex = pixel2texpair(iw, ih, w, h);

        let (ray_org, ray_dir) = del_raycast_core::cam_pbrt::cast_ray(
            (iw, ih),
            (0., 0.),
            (w, h),
            fov,
            transform_cam_lcl2glbl,
        );

        let mut halton = HaltonSampler::new();
        let nsamples = 64;

        let mut result = [0.; 3];
        for i in 0..nsamples {
            todo!()
        }

        pix.0 = result;
    };

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix));

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/dof.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&img, w, h);
}
