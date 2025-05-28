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

fn spherical_harmonics(l: i32, m: i32, theta: f32, phi: f32) -> f32 {
    todo!()
}

fn main() {
    use image::Pixel;

    let w = 512;
    let h = 512;
    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    // center, radius
    let sphere = ([0f32, 0f32, 0f32], 1f32);

    // TODO: use y-up
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
            let spherical_coord = xyz2spherical_local(&hit_pos, &sphere);

            result = [spherical_coord[2].sin(); 3];
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

#[test]
fn test() {}
