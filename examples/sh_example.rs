use std::f32::consts::PI;

use del_geo_core::vec3::{self, Vec3};
use rayon::prelude::*;
use rs_sampler::cam;

type Rgb = image::Rgb<f32>;

/// y-up, z-forward
fn xyz2spherical_local(pos: &[f32; 3], sphere: &([f32; 3], f32)) -> [f32; 3] {
    let v = vec3::sub(pos, &sphere.0);
    let r = v.norm();
    assert!(r > 0f32);
    let u = v.normalize();
    let theta = u[1].acos();
    let phi = u[0].atan2(u[2]);
    [r, theta, phi]
}

/// evaluate an Associated Legendre Polynomial P(l,m) at x, using three properties
#[allow(non_snake_case)]
fn P(l: i32, m: i32, x: f32) -> f32 {
    assert!(l >= m.abs());
    let mut pmm = 1f32;
    // evaluate  P(m,m) from P(0,0)
    if m > 0 {
        let sqrtfactor = ((1f32 - x) * (1f32 + x)).sqrt();
        let mut fact = 1f32;
        for _ in 1..m + 1 {
            pmm *= -fact * sqrtfactor;
            fact += 2f32;
        }
    }
    if l == m {
        return pmm;
    }

    let mut pmm1 = x * (2f32 * m as f32 + 1f32) * pmm;
    if l == m + 1 {
        return pmm1;
    }

    let mut pll = 0f32;
    for ll in m + 2..l + 1 {
        pll = (x * (2 * ll - 1) as f32 * pmm1 - (ll + m - 1) as f32 * pmm) / (ll - m) as f32;
        pmm = pmm1;
        pmm1 = pll;
    }

    pll
}

fn factorial(x: i32) -> i32 {
    if x == 0 {
        return 1;
    }
    let mut acc = 1;
    let mut n = x;
    while n > 0 {
        acc *= n;
        n -= 1;
    }

    acc
}

#[allow(non_snake_case)]
fn K(l: i32, m: i32) -> f32 {
    let mabs = m.abs();
    let fac0 = (2f32 * l as f32 + 1f32) / (4f32 * PI);
    let fa1 = factorial(l - mabs) as f32;
    let fa2 = factorial(l + mabs) as f32;
    let fac1 = fa1 / fa2;
    let res = fac0 * fac1;
    res.sqrt()
}

/// real part of spherical harmonics
fn sh_real(l: i32, m: i32, theta: f32, phi: f32) -> f32 {
    if m == 0 {
        return K(l, m) * P(l, m, theta.cos());
    }

    let sqrt2 = 2f32.sqrt();

    if m > 0 {
        let k = K(l, m);
        let p = P(l, m, theta.cos());
        return sqrt2 * k * (m as f32 * phi) * p;
        // return sqrt2 * K(l, m) * (m as f32 * phi) * P(l, m, theta.cos());
    }

    let m = -m;
    sqrt2 * K(l, m) * (m as f32 * phi) * P(l, m, theta.cos())
}

fn draw_legendre_poly() {
    use image::Pixel;

    let w = 512;
    let h = 512;
    let mut img = vec![*Rgb::from_slice(&[0f32; 3]); w * h];

    let deg = 5;
    let ord = 0;
    for l in 0..deg {
        for i in 0..w {
            let x = i as f32 / w as f32;
            let y = P(l, ord, 2f32 * x - 1f32);
            let iy = ((1. - (y + 1f32) / 2f32) * h as f32 - 0.5f32) as usize;
            let ix = (x * w as f32 - 0.5f32) as usize;
            let ipix = iy * w + ix;
            img[ipix].0 = [1f32; 3];
        }
    }

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/sh_example_legendre.hdr").unwrap();
    let _ = HdrEncoder::new(file_ms).encode(&img, w, h);
}

fn draw_sh() {
    use image::Pixel;

    let w = 512;
    let h = 512;
    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    // center, radius
    let sphere = ([0f32, 0f32, 0f32], 1f32);

    let campos = [0., 0., 3.];
    let view = [0., 0., -1.];
    let mut v2w = cam::matrix_v2w(&view).1;
    // concat translation
    v2w[0 + 3 * 4] = campos[0];
    v2w[1 + 3 * 4] = campos[1];
    v2w[2 + 3 * 4] = campos[2];

    let task = |i_pix: usize, pix: &mut Rgb| {
        let iw = i_pix % w;
        let ih = i_pix / w;
        let mut result = [0f32; 3];

        let (ray_org, ray_dir) = cam::gen_ray((iw, ih), (0., 0.), (w, h), 60., &v2w);

        let hit_res = del_geo_nalgebra::sphere::intersection_ray(
            &nalgebra::Vector3::<f32>::from(sphere.0),
            sphere.1,
            &nalgebra::Vector3::<f32>::from(ray_org),
            &nalgebra::Vector3::<f32>::from(ray_dir),
        );

        if let Some(t) = hit_res {
            let hit_pos = vec3::axpy::<f32>(t, &ray_dir, &ray_org);
            let coord = xyz2spherical_local(&hit_pos, &sphere);
            let theta = coord[1];
            let phi = coord[2];

            let v = sh_real(1, 1, theta, phi);

            if v > 0f32 {
                result = [v, 0f32, 0f32];
            } else {
                result = [0f32, -v, 0f32];
            }
        }

        pix.0 = result;
    };

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| task(i_pix, pix));

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/sh_example.hdr").unwrap();
    let _ = HdrEncoder::new(file_ms).encode(&img, w, h);
}

fn main() {
    draw_sh();
    draw_legendre_poly();
}

#[test]
fn test() {}
