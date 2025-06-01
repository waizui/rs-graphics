use del_geo_core::vec3::{self, Vec3};
use rayon::prelude::*;
use rs_sampler::cam;
use std::f32::consts::PI;

type Rgb = image::Rgb<f32>;

struct Sphere {
    cnt: [f32; 3],
    r: f32,
}

/// y-up, z-forward
fn xyz2spherical(xyz: &[f32; 3]) -> [f32; 3] {
    let r = xyz.norm();
    assert!(r > 0f32);
    let u = xyz.normalize();
    let theta = u[1].acos();
    let phi = u[0].atan2(u[2]);
    [r, theta, phi]
}

/// y-up, z-forward
/// return local spherical coordinates of sphere
fn worldpos2sphere(pos: &[f32; 3], sphere: &Sphere) -> [f32; 3] {
    let v = vec3::sub(pos, &sphere.cnt);
    xyz2spherical(&v)
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
    let fac1 = factorial(l - mabs) as f32 / factorial(l + mabs) as f32;
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
        return sqrt2 * K(l, m) * (m as f32 * phi).cos() * P(l, m, theta.cos());
    }

    let m = -m;
    sqrt2 * K(l, m) * (m as f32 * phi).sin() * P(l, m, theta.cos())
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

/// return (sphere,l,m)
fn gen_spheres(l: i32) -> Vec<(Sphere, (i32, i32))> {
    let mut spheres: Vec<(Sphere, (i32, i32))> = Vec::new();
    let ext = 3;
    for il in 0..l {
        for im in -il..il + 1 {
            let x = (im * ext) as f32;
            let y = ((l - il) * ext) as f32 - (l * ext) as f32 / 2f32;
            let s = Sphere {
                cnt: [x, y, 0f32],
                r: 1f32,
            };

            spheres.push((s, (il, im)));
        }
    }
    spheres
}

fn draw_sh() {
    use image::Pixel;

    let w = 512;
    let h = 512;
    let mut img = vec![*Rgb::from_slice(&[0.5; 3]); w * h];

    let l = 4; //sh order
    let spheres = gen_spheres(l);

    let campos = [0., 0., 18.];
    let view = [0., 0., -1.];
    let mut v2w = cam::matrix_v2w(&view).1;
    // concat translation
    v2w[0 + 3 * 4] = campos[0];
    v2w[1 + 3 * 4] = campos[1];
    v2w[2 + 3 * 4] = campos[2];

    let task = |i_pix: usize, pix: &mut Rgb| {
        let iw = i_pix % w;
        let ih = i_pix / w;
        let (ray_org, ray_dir) = cam::gen_ray((iw, ih), (0., 0.), (w, h), 60., &v2w);

        let (t, i_sphere) = {
            let mut t = f32::INFINITY;
            let mut i_sphere = 0;
            for (i, sphere) in spheres.iter().enumerate() {
                if let Some(t_hit) = del_geo_nalgebra::sphere::intersection_ray(
                    &nalgebra::Vector3::<f32>::from(sphere.0.cnt),
                    sphere.0.r,
                    &nalgebra::Vector3::<f32>::from(ray_org),
                    &nalgebra::Vector3::<f32>::from(ray_dir),
                ) {
                    if t_hit < t {
                        t = t_hit;
                        i_sphere = i;
                    }
                }
            }

            (t, i_sphere)
        };

        if t == f32::INFINITY {
            return;
        }

        let sphere = &spheres[i_sphere];
        let hit_pos = vec3::axpy::<f32>(t, &ray_dir, &ray_org);
        let coord = worldpos2sphere(&hit_pos, &sphere.0);
        let theta = coord[1];
        let phi = coord[2];

        let v = sh_real(sphere.1 .0, sphere.1 .1, theta, phi);

        if v > 0f32 {
            pix.0 = [0f32, v, 0f32];
        } else {
            pix.0 = [-v, 0f32, 0f32];
        }
    };

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| task(i_pix, pix));

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/sh_example_spheres.hdr").unwrap();
    let _ = HdrEncoder::new(file_ms).encode(&img, w, h);
}

#[derive(Clone, Debug)]
struct SHSample {
    sph: [f32; 3],
    xyz: [f32; 3],
    coeff: Vec<f32>,
}

/// nsamples: specify how many samples will be generated
/// return random sh samples  across sphere surface
fn gen_sh_samples(nsamples: usize, l: i32) -> Vec<SHSample> {
    use rs_sampler::haltonsampler::radical_inverse;
    use rs_sampler::sampling::*;

    let mut samples = vec![
        SHSample {
            sph: [0f32; 3],
            xyz: [0f32; 3],
            coeff: Vec::new()
        };
        nsamples
    ];

    let task = |isample: usize, sample: &mut SHSample| {
        // generate random samples
        let rx: f32 = radical_inverse(isample, 2);
        let ry: f32 = radical_inverse(isample, 3);
        let xyz = square2unitsphere(&[rx, ry]);
        let spherial = xyz2spherical(&xyz);
        let theta = spherial[1];
        let phi = spherial[2];

        sample.xyz = xyz;
        sample.sph = [1f32, theta, phi];

        for il in 0..l + 1 {
            for im in -il..il + 1 {
                let sh = sh_real(il, im, theta, phi);
                sample.coeff.push(sh);
            }
        }
    };

    samples
        .par_iter_mut()
        .enumerate()
        .for_each(|(isample, sample)| task(isample, sample));

    samples
}

fn project_sh() {
    use image::Pixel;

    let nsamples = 128;
    let l = 3;
    let sh_samples = gen_sh_samples(nsamples, l);

    let w = 512;
    let h = 512;
    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    let sphere = Sphere {
        cnt: [0f32; 3],
        r: 1f32,
    };

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
        let (ray_org, ray_dir) = cam::gen_ray((iw, ih), (0., 0.), (w, h), 60., &v2w);

        if let Some(t) = del_geo_nalgebra::sphere::intersection_ray(
            &nalgebra::Vector3::<f32>::from(sphere.cnt),
            sphere.r,
            &nalgebra::Vector3::<f32>::from(ray_org),
            &nalgebra::Vector3::<f32>::from(ray_dir),
        ) {
            let hit_pos = vec3::axpy::<f32>(t, &ray_dir, &ray_org);
            let coord = worldpos2sphere(&hit_pos, &sphere);
            let theta = coord[1];
            let phi = coord[2];

            pix.0 = [phi.cos().abs(); 3];
            for sh_sample in &sh_samples {
                let sh_theta = sh_sample.sph[1];
                if (sh_theta.abs() - theta.abs()).abs() > 1e-2 {
                    continue;
                }

                let sh_phi = sh_sample.sph[2];
                if (sh_phi.abs() - phi.abs()).abs() < 1e-2 {
                    pix.0 = [1f32, 0f32, 0f32];
                }
            }
        };
    };

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| task(i_pix, pix));

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/sh_example_projection.hdr").unwrap();
    let _ = HdrEncoder::new(file_ms).encode(&img, w, h);
}

fn main() {
    draw_legendre_poly();
    draw_sh();
    project_sh();
}

#[test]
fn test() {}
