use del_raycast_core::io_pfm::PFM;
use rs_sampler::{haltonsampler::HaltonSampler, sampler::Sampler};

type Image = PFM;
type Real = f32;

fn calc_gray_scale(path: &str) -> Image {
    let pfm = PFM::read_from(path).unwrap();

    if pfm.channels < 3 {
        panic!("pfm channel invalid");
    }

    let mut gray_scale = PFM {
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
            gray_scale.data[h * pfm.w + w] = gray;
        }
    }

    gray_scale
}

fn roundi(f: Real) -> i32 {
    (f + 0.5) as i32
}

fn tex2pixel(t: Real, num_p: usize) -> usize {
    roundi(t * (num_p as Real) - 0.5) as usize
}

/// x,y: uniformal variables
/// returns sampling position of x,y
fn calc_sample_coord(x: Real, y: Real, img: &Image) -> (usize, usize) {
    let w = img.w;
    let h = img.h;
    let mut sum = 0.;

    let mut r_avg = 0.;
    for i_h in 0..h {
        let i_w = tex2pixel(x, w); // at given y
        let gray = img.data[i_w * w + i_h];
        r_avg += gray;
    }

    r_avg /= w as Real;

    let mut i_x = 0; // sampling pixel coord of variable y
    for i_w in 0..w {
        let i_h = tex2pixel(y, h); // at given y
        let gray = img.data[i_h * w + i_w];
        sum += gray / r_avg; // TODO: r_avg ==0
        if sum >= y {
            i_x = i_w;
            break;
        }
    }

    let mut c_avg = 0.;
    for i_w in 0..w {
        let i_h = tex2pixel(y, h);
        let gray = img.data[i_w * w + i_h];
        c_avg += gray;
    }
    c_avg /= w as Real;

    sum = 0.;
    let mut i_y = 0;
    // col sum
    for i_h in 0..h {
        let i_w = tex2pixel(x, w); // at given x
        let gray = img.data[i_h * w + i_w];
        sum += gray / c_avg;
        if sum >= x {
            i_y = i_h;
            break;
        }
    }

    (i_x, i_y)
}

fn sample_light(img: &Image) -> Image {
    let w = img.w;
    let h = img.h;

    let mut res = PFM {
        data: vec![0.; w * h],
        w,
        h,
        channels: 1,
        little_endian: img.little_endian,
    };

    let mut halton = HaltonSampler::new();

    for i_h in 0..h {
        for i_w in 0..w {
            Sampler::<Real>::set_i(&mut halton, i_h * w + i_w);
            Sampler::<Real>::set_dim(&mut halton, 0);
            let xy: [Real; 2] = halton.get2d();
            let st = calc_sample_coord(xy[0], xy[1], img);
            res.data[st.0 * h + st.1] = 1.;
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
    let img = calc_gray_scale("asset/envmap.pfm");
    let env = sample_light(&img);
    let _ = env.save("target/envmap_gray.pfm");
}
