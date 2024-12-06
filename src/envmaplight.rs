type Real = f32;
type Rgb = image::Rgb<Real>;

use image::Pixel;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;

pub fn calc_grayscale(img: &Vec<Rgb>, w: usize, h: usize) -> Vec<Rgb> {
    let iter = |i_pix: usize, pix: &mut Rgb| {
        let c = &img[i_pix];
        let gray = (0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]).clamp(0., 1.);
        pix.0 = [gray; 3];
    };

    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix));

    img
}

#[test]
fn test_grayscale() {
    use crate::pfm::PFM;
    use itertools::Itertools;
    let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    let rgbdata = pfm
        .data
        .chunks(pfm.channels)
        .map(|chunk| *Rgb::from_slice(&[chunk[0], chunk[0], chunk[0]]))
        .collect_vec();

    let grayscale = calc_grayscale(&rgbdata, pfm.w, pfm.h);
    let file_ms = std::fs::File::create("target/04_env_light_grayscale.hdr").unwrap();
    use image::codecs::hdr::HdrEncoder;
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&grayscale, pfm.w, pfm.h);
}
