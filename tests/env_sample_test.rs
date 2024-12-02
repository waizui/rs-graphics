use del_raycast_core::io_pfm::PFM;
use rs_sampler::{haltonsampler::HaltonSampler, sampler::Sampler};

type Image = PFM;
type Real = f32;
const PI: Real = std::f32::consts::PI;

fn roundi(f: Real) -> i32 {
    (f + 0.5) as i32
}

fn tex2pixel(t: Real, num_p: usize) -> usize {
    roundi(t * (num_p as Real) - 0.5) as usize
}

fn pixel2tex(x: usize, y: usize, img: &Image) -> [Real; 2] {
    let u = x as Real / img.w as Real;
    let v = y as Real / img.h as Real;
    [u, v]
}

fn unitsphere2envmap(d: &[f32; 3]) -> [Real; 2] {
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
    let v = v.copysign(-d[1]);
    [u * 0.5 + 0.5, v * 0.5 + 0.5]
}

fn envmap2unitsphere(u: Real, v: Real) -> [Real; 3] {
    let r = u - 0.5;
    let phi = (PI / 4.) * (v - 0.5) / r;

    let x = phi.cos() * r * (2. - r * r).sqrt();
    let y = phi.sin() * r * (2. - r * r).sqrt();
    let z = 1. - r * r;

    [x, y, z]
}

fn get_color(u: Real, v: Real, img: &Image) -> [Real; 3] {
    let channel = img.channels;
    let i_w = tex2pixel(u, img.w);
    let i_h = tex2pixel(v, img.w);

    // let r = img.data[530360];
    // let g = img.data[530361];
    // let b = img.data[530362];
    //
    // dbg!(r,g,b);

    let i = i_h * img.w + i_w;
    if channel == 3 {
        let r = img.data[i];
        let g = img.data[i + 1];
        let b = img.data[i + 2];
        [r, g, b]
    } else if channel == 1 {
        let g = img.data[i];
        [g, g, g]
    } else {
        panic!("Invalid channel");
    }
}

fn set_color(col: &[Real; 3], u: Real, v: Real, img: &mut Image) {
    let channel = img.channels;
    let i_w = tex2pixel(u, img.w);
    let i_h = tex2pixel(v, img.w);

    let i = i_h * img.w + i_w;
    if channel == 3 {
        img.data[i] = col[0];
        img.data[i + 1] = col[1];
        img.data[i + 2] = col[2];
    } else if channel == 1 {
        img.data[i] = col[0];
    } else {
        panic!("Invalid channel");
    }
}

fn calc_grayscale(pfm: &Image) -> Image {
    if pfm.channels < 3 {
        panic!("pfm channel invalid");
    }

    let mut grayscale = PFM {
        data: vec![0.; pfm.w * pfm.h],
        w: pfm.w,
        h: pfm.h,
        channels: 1,
        little_endian: pfm.little_endian,
    };

    let channels = pfm.channels;
    for h in 0..pfm.h {
        for w in 0..pfm.w {
            let p_i = channels * (h * pfm.w + w);
            let r = pfm.data[p_i];
            let g = pfm.data[p_i + 1];
            let b = pfm.data[p_i + 2];
            let gray = (0.2126 * r + 0.7152 * g + 0.0722 * b).clamp(0., 1.);
            grayscale.data[h * pfm.w + w] = gray;
        }
    }

    grayscale
}

fn calc_integral_over_grayscale(img: &Image) -> Real {
    let w = img.w;
    let h = img.h;

    let mut sum = 0.;

    for i_h in 0..h {
        for i_w in 0..w {
            let gray = img.data[i_h * w + i_w];
            sum += gray;
        }
    }

    sum /= (w * h) as Real;
    sum
}

fn binary_find_inv_cdf() -> Real {
    1.
}

/// x,y: uniformal variables
/// returns sampling position of x,y and pdf
fn calc_sample_coord(x: Real, y: Real, img: &Image, itgr: Real) -> ([usize; 2], Real) {
    let w = img.w;
    let h = img.h;
    let mut sum = 0.;

    let mut r_avg = 0.;
    let i_h_fix = tex2pixel(y, h);
    for i_w in 0..w {
        let gray = img.data[i_h_fix * w + i_w];
        r_avg += gray;
    }
    r_avg /= w as Real;

    let mut i_x = 0; // sampling pixel coord of variable y
    for i_w in 0..w {
        let gray = img.data[i_h_fix * w + i_w];
        sum += gray / (r_avg * w as Real); // TODO: r_avg ==0
        if sum >= x {
            i_x = i_w;
            break;
        }
    }

    let mut c_avg = 0.;
    let i_w_fix = tex2pixel(x, w); // at given y
    for i_h in 0..h {
        let gray = img.data[i_h * w + i_w_fix];
        c_avg += gray;
    }
    c_avg /= h as Real;

    sum = 0.;
    let mut i_y = 0;
    for i_h in 0..h {
        let gray = img.data[i_h * w + i_w_fix];
        sum += gray / (c_avg * h as Real);
        if sum >= y {
            i_y = i_h;
            break;
        }
    }

    let pdf = img.data[i_h_fix * w + i_w_fix] / itgr;

    ([i_x, i_y], pdf)
}

fn sample_light(grayscale: &Image, envmap: &Image) -> Image {
    let w = grayscale.w;
    let h = grayscale.h;

    let mut res = PFM {
        data: vec![0.; w * h],
        w,
        h,
        channels: 3,
        little_endian: grayscale.little_endian,
    };

    let mut halton = HaltonSampler::new();
    let itgr = calc_integral_over_grayscale(grayscale);

    let nsamples = 8;
    for i_w in 0..64 {
        for i_h in 0..64 {
            // screen coord to world
            let tex = pixel2tex(i_w, i_h, grayscale);
            let nrm = envmap2unitsphere(tex[0], tex[1]);

            for i in 0..nsamples {
                Sampler::<Real>::set_i(&mut halton, i + 1);
                Sampler::<Real>::set_dim(&mut halton, 0);
                let xy: [Real; 2] = halton.get2d();
                // sampling position
                let coord_pfd = calc_sample_coord(xy[0], xy[1], grayscale, itgr);
                let st = coord_pfd.0;
                let p2t = pixel2tex(st[0], st[1], grayscale);
                // sampling ray direction
                let raydir = envmap2unitsphere(p2t[0], p2t[1]);

                let costheta = del_geo_core::vec3::dot(&raydir, &nrm);
                if costheta > 0. {
                    let mut color = get_color(p2t[0], p2t[1], envmap);
                    let pdf = coord_pfd.1;
                    set_color(&color, tex[0], tex[1], &mut res);
                }
                // costheta
            }
        }
    }

    res
}

fn calc_cdf_inv_lookup(img: &Image) -> [Image; 2] {
    todo!()
}

fn sample_light_lookup(img: &Image, tbl: &[Image; 2]) -> Image {
    todo!()
}

#[test]
fn test_env_sample() {
    let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    let grayscale = calc_grayscale(&pfm);
    let env = sample_light(&grayscale, &pfm);
    let _ = env.save("target/envmap_gray.pfm");
}
