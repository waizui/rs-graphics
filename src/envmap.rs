type Real = f32;
type Rgb = image::Rgb<Real>;

use image::Pixel;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;

pub fn roundi(f: Real) -> i32 {
    (f + 0.5) as i32
}

/// t:texel, ext: extent
/// texel coordinates([0,1]^2) to pixel coordinates ([w,h])
pub fn tex2pixel(t: Real, ext: usize) -> usize {
    roundi(t * ((ext - 1) as Real) - 0.5) as usize
}

/// p: index of pixel, ext: extent
/// pixel coordinates ([w,h]) to texel coordinates([0,1]^2)
pub fn pixel2tex(p: usize, ext: usize) -> Real {
    if ext != 1 {
        return p as Real / (ext - 1) as Real;
    }

    0.
}

/// pixel coordinates ([w,h]) to texel coordinates ([0,1]^2)
pub fn pixel2texpair(x: usize, y: usize, w: usize, h: usize) -> [Real; 2] {
    let u = pixel2tex(x, w);
    let v = pixel2tex(y, h);
    [u, v]
}

pub fn calc_grayscale(img: &[Rgb], w: usize, h: usize) -> Vec<Rgb> {
    let iter = |i_pix: usize, pix: &mut Rgb| {
        let c = img[i_pix];
        let gray = 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];
        pix.0 = [gray; 3];
    };

    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix));

    img
}

/// return columns average map and rows average map
pub fn calc_avg_col_row(grayscale: &[Rgb], w: usize, h: usize) -> [Vec<Rgb>; 2] {
    let iter = |i_pix: usize, pix: &mut Rgb, is_row_avg: bool| {
        let ext = if is_row_avg { w } else { h };

        let mut avg = 0.;
        for i in 0..ext {
            let gray = if is_row_avg {
                // elements average of i_pix-th row
                grayscale[i_pix * w + i][0]
            } else {
                // elements average of i_pix-th column
                grayscale[i * w + i_pix][0]
            };

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
        .for_each(|(i_pix, pix)| iter(i_pix, pix, true));

    [col_avg_map, row_avg_map]
}

/// return marginal probabilty map of y and conditional probabilty map of p(x|y)
pub fn calc_inverse_cdf_map(
    grayscale: &[Rgb],
    itgr: Real,
    w: usize,
    h: usize,
) -> (Vec<Rgb>, Vec<Rgb>) {
    let avgs = calc_avg_col_row(grayscale, w, h);
    let avg = &avgs[1];

    let build_marginal_map = |i_pix: usize, pix: &mut Rgb| {
        let cur_h = i_pix;
        let mut sum = 0.;
        let mut i_y = 0;
        let v = pixel2tex(cur_h, h);

        for (i_h, item) in avg.iter().enumerate().take(h) {
            // row_avg(h)/itgr
            sum += item[0] / (h as Real * itgr);
            if sum >= v {
                i_y = i_h;
                break;
            }
        }

        let invy = pixel2tex(i_y, h);
        pix.0 = [invy; 3];
    };

    let mut marginal_map = vec![*Rgb::from_slice(&[0.; 3]); h];
    marginal_map
        .par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| build_marginal_map(i_pix, pix));

    let mut conditinal_map = vec![*Rgb::from_slice(&[0.; 3]); w * h];
    let build_conditional_map = |i_pix: usize, pix: &mut Rgb| {
        let mut sum = 0.;

        let cur_w = i_pix % w;
        let cur_h = i_pix / w;

        let mut i_x = 0;
        let u = pixel2tex(cur_w, w);
        for i_w in 0..w {
            // p(w|h) = p(w,h)/p(h) = gray(w,h) / row_avg(w)
            let avg = avg[cur_h][0];
            sum += grayscale[cur_h * w + i_w][0] / (avg * w as Real);
            if sum >= u {
                i_x = i_w;
                break;
            }
        }

        let invx = pixel2tex(i_x, w);
        pix.0 = [invx; 3];
    };

    conditinal_map
        .par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| build_conditional_map(i_pix, pix));

    (marginal_map, conditinal_map)
}

pub fn calc_integral_over_grayscale(grayscale: &[Rgb], w: usize, h: usize) -> Real {
    let mut sum = 0.;

    for i_h in 0..h {
        for i_w in 0..w {
            let c = grayscale[i_w + i_h * w][0];
            sum += c;
        }
    }

    sum /= (w * h) as Real;
    sum
}

/// d: xyz at sphere surface
pub fn unitsphere2envmap(d: &[f32; 3]) -> [Real; 2] {
    let x = d[0].abs();
    let y = d[1].abs();
    let z = d[2].abs();
    let r = (1. - z).sqrt();
    let phi = y.atan2(x);
    let phi = phi * std::f32::consts::FRAC_2_PI;
    let v = phi * r;
    let u = r - v;
    let (u, v) = if d[2] < 0. { (1. - v, 1. - u) } else { (u, v) };
    let u = u.copysign(d[0]);
    let v = v.copysign(d[1]);
    [u * 0.5 + 0.5, v * 0.5 + 0.5]
}

// https://github.com/mmp/pbrt-v4/blob/1ae72cfa7344e79a7815a21ed3da746cdccee59b/src/pbrt/util/math.cpp#L292
pub fn envmap2unitsphere(p: &[Real; 2]) -> [Real; 3] {
    use std::f32::consts::PI;

    // map [0,1] to [-1,1]
    let u = 2. * p[0] - 1.;
    let v = 2. * p[1] - 1.;
    let up = u.abs();
    let vp = v.abs();
    let sd = 1. - (up + vp);
    let d = sd.abs();
    let r = 1. - d;
    let phi = (if r == 0. { 1. } else { (vp - up) / r + 1. }) * PI / 4.;
    let z = (1. - r * r).copysign(sd);
    let cosphi = phi.cos().copysign(u);
    let sinphi = phi.sin().copysign(v);

    let x = cosphi * r * (2. - r * r).sqrt();
    let y = sinphi * r * (2. - r * r).sqrt();

    [x, y, z]
}
#[test]
fn test_sphere_mapping() {
    let mut dir = [0.; 3];

    for i in 0..100 {
        for j in 0..3 {
            dir[j] = (i + 1) as Real;
        }

        del_geo_core::vec3::normalize(&mut dir);

        let uv = unitsphere2envmap(&dir);
        let udir = envmap2unitsphere(&uv);

        assert!((dir[0] - udir[0]).abs() < 0.001);
        assert!((dir[1] - udir[1]).abs() < 0.001);
        assert!((dir[2] - udir[2]).abs() < 0.001);
    }
}

#[test]
fn test_grayscale() {
    use crate::pfm::PFM;
    use itertools::Itertools;
    let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    let rgbdata = pfm
        .data
        .chunks(pfm.channels)
        .map(|chunk| *Rgb::from_slice(&[chunk[0], chunk[1], chunk[2]]))
        .collect_vec();

    let grayscale = calc_grayscale(&rgbdata, pfm.w, pfm.h);
    let file_ms = std::fs::File::create("target/envmap_light_grayscale.hdr").unwrap();
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
        .map(|chunk| *Rgb::from_slice(&[chunk[0], chunk[1], chunk[2]]))
        .collect_vec();
    let grayscale = calc_grayscale(&rgbdata, pfm.w, pfm.h);

    let avg = calc_avg_col_row(&grayscale, pfm.w, pfm.h);

    use image::codecs::hdr::HdrEncoder;

    let file_ms = std::fs::File::create("target/envmap_light_row_avg.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&avg[1], 1, pfm.h);

    let file_ms = std::fs::File::create("target/envmap_light_col_avg.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&avg[0], pfm.w, 1);
}

#[test]
fn test_inverse_cdf() {
    use crate::haltonsampler::HaltonSampler;
    use crate::pfm::PFM;
    use crate::sampler::Sampler;
    use itertools::Itertools;
    use std::time::Instant;

    let mut sw = Instant::now();

    let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    let w = pfm.w;
    let h = pfm.h;
    let mut rgbdata = pfm
        .data
        .chunks(pfm.channels)
        .map(|chunk| *Rgb::from_slice(&[chunk[0], chunk[1], chunk[2]]))
        .collect_vec();

    let grayscale = calc_grayscale(&rgbdata, pfm.w, pfm.h);
    let itgr = calc_integral_over_grayscale(&grayscale, pfm.w, pfm.h);
    let (marginal_map, conditional_map) = calc_inverse_cdf_map(&grayscale, itgr, w, h);

    println!("Build inverse cdf map:{} ms", sw.elapsed().as_millis());

    sw = Instant::now();

    let mut halton = HaltonSampler::new();
    for i in 0..1024 * 1024 {
        Sampler::<Real>::set_i(&mut halton, i);
        Sampler::<Real>::set_dim(&mut halton, 0);
        let random: [Real; 2] = Sampler::get2d(&mut halton);
        let r_w = tex2pixel(random[0], w);
        let r_h = tex2pixel(random[1], h);

        let sampley = marginal_map[r_h][0];
        let samplex = conditional_map[r_w + tex2pixel(sampley, h) * w][0];

        let p_w = tex2pixel(samplex, w);
        let p_h = tex2pixel(sampley, h);
        rgbdata[p_h * w + p_w] = *Rgb::from_slice(&[1., 0., 0.]);
    }

    println!("Generate Samples:{} ms", sw.elapsed().as_millis());

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/envmap_light_invcdf.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&rgbdata, w, h);
}
