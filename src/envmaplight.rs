type Real = f32;
type Rgb = image::Rgb<Real>;

use image::Pixel;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;

pub fn calc_grayscale(img: &Vec<Rgb>, w: usize, h: usize) -> Vec<Rgb> {
    let iter = |i_pix: usize, pix: &mut Rgb| {
        let c = img[i_pix];
        let gray = (0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]).clamp(0., 1.);
        pix.0 = [gray; 3];
    };

    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix));

    img
}

pub fn calc_avg_col_row(grayscale: &Vec<Rgb>, w: usize, h: usize) -> [Vec<Rgb>; 2] {
    let iter = |i_pix: usize, pix: &mut Rgb, is_row: bool| {
        let (stride, ext) = if is_row {
            // calculate column averages for a row
            // march w element to get same column element in next row
            (w, h)
        } else {
            // get next element in same row
            (1, w)
        };

        let mut avg = 0.;
        for i in 0..ext {
            let gray = grayscale[i_pix + i * stride][0];
            avg += gray;
        }
        avg /= ext as Real;
        pix.0 = [avg; 3];
    };

    let mut col_avg_map = vec![*Rgb::from_slice(&[0.; 3]); w];
    let mut row_avg_map = vec![*Rgb::from_slice(&[0.; 3]); h];

    col_avg_map
        .par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix, false));

    row_avg_map
        .par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix, false));

    [col_avg_map, row_avg_map]
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

#[test]
fn test_calc_col_row_avg() {
    use crate::pfm::PFM;
    use itertools::Itertools;
    let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    let rgbdata = pfm
        .data
        .chunks(pfm.channels)
        .map(|chunk| *Rgb::from_slice(&[chunk[0], chunk[0], chunk[0]]))
        .collect_vec();
    let grayscale = calc_grayscale(&rgbdata, pfm.w, pfm.h);

    let avg = calc_avg_col_row(&grayscale, pfm.w, pfm.h);

    use image::codecs::hdr::HdrEncoder;

    let file_ms = std::fs::File::create("target/04_env_light_col_avg.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&avg[0], 1, pfm.h);

    let file_ms = std::fs::File::create("target/04_env_light_row_avg.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&avg[1], pfm.w, 1);
}
